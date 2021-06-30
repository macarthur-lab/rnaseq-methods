import argparse
import hail as hl  # used for hadoop file utils
import hailtop.batch as hb
import logging
import os

from sample_metadata.rnaseq_metadata_utils import get_rnaseq_metadata_joined_with_paths_df

hl.init(log="/dev/null")

logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s', level=logging.INFO)
logger = logging.getLogger(__name__)


GCLOUD_PROJECT = "seqr-project"
GCLOUD_USER_ACCOUNT = "weisburd@broadinstitute.org"
GCLOUD_CREDENTIALS_LOCATION = "gs://weisburd-misc/creds"

def main():
    rnaseq_sample_metadata_df = get_rnaseq_metadata_joined_with_paths_df()

    p = argparse.ArgumentParser()
    grp = p.add_mutually_exclusive_group(required=True)
    grp.add_argument("--local", action="store_true", help="Batch: run locally")
    grp.add_argument("--cluster", action="store_true", help="Batch: submit to cluster")
    p.add_argument("--batch-billing-project", default="tgg-rare-disease", help="Batch: billing project. Required if submitting to cluster.")
    p.add_argument("--batch-job-name", help="Batch: (optional) job name")

    p.add_argument("-f", "--force", action="store_true", help="Recompute and overwrite cached or previously computed data")
    grp = p.add_mutually_exclusive_group(required=True)
    grp.add_argument("-b", "--rnaseq-batch-name", nargs="*", help="RNA-seq batch names to process", choices=set(rnaseq_sample_metadata_df['star_pipeline_batch']))
    grp.add_argument("-s", "--rnaseq-sample-id", nargs="*", help="RNA-seq sample IDs to process", choices=set(rnaseq_sample_metadata_df['sample_id']))
    args = p.parse_args()

    #logger.info("\n".join(df.columns))

    if args.rnaseq_batch_name:
        batch_names = args.rnaseq_batch_name
        sample_ids = rnaseq_sample_metadata_df[rnaseq_sample_metadata_df['star_pipeline_batch'].isin(batch_names)].sample_id
    elif args.rnaseq_sample_id:
        sample_ids = args.rnaseq_sample_id

    logger.info(f"Processing {len(sample_ids)} sample ids: {', '.join(sample_ids)}")

    # see https://hail.zulipchat.com/#narrow/stream/223457-Batch-support/topic/auth.20as.20user.20account for more details
    if args.local:
        backend = hb.LocalBackend(gsa_key_file=os.path.expanduser("~/.config/gcloud/misc-270914-cb9992ec9b25.json"))
    else:
        backend = hb.ServiceBackend(args.batch_billing_project)

    b = hb.Batch(backend=backend, name=args.batch_job_name)

    # define workflow inputs
    if args.local:
        genes_gtf = b.read_input("gencode.v26.annotation.gff3", extension=".gff3")
    else:
        genes_gtf = b.read_input("gs://tgg-rnaseq/ref/gencode.v26.annotation.GRCh38.gff3", extension=".gff3")

    # define parallel execution for samples
    for sample_id in sample_ids:
        metadata_row = rnaseq_sample_metadata_df.loc[sample_id]
        batch_name = metadata_row['star_pipeline_batch']

        # set job inputs & outputs
        input_read_data = b.read_input_group(
            bam=metadata_row['star_bam'],
            bai=metadata_row['star_bai'],
        )

        output_dir = f"gs://tgg-rnaseq/{batch_name}/majiq_build/"
        output_file_path = os.path.join(output_dir, f"majiq_build_{sample_id}.tar.gz")

        # check if output file already exists
        if hl.hadoop_is_file(output_file_path) and not args.force:
            logger.info(f"{sample_id} output file already exists: {output_file_path}. Skipping...")
            continue

        file_stats = hl.hadoop_stat(metadata_row['star_bam'])
        bam_size = int(round(file_stats['size_bytes']/10.**9))

        # define majiq build commands for this sample
        j = b.new_job(name=args.batch_job_name)
        j.image("weisburd/majiq:latest")
        j.storage(f'{bam_size*3}Gi')
        j.cpu(1)  # default: 1
        j.memory("15G")  # default: 3.75G
        logger.info(f'Requesting: {j._storage or "default"} storage, {j._cpu or "default"} CPU, {j._memory or "default"} memory')

        # switch to user account
        j.command(f"gcloud auth activate-service-account --key-file /gsa-key/key.json")
        j.command(f"gsutil -m cp -r {GCLOUD_CREDENTIALS_LOCATION}/.config /tmp/")
        j.command(f"rm -rf ~/.config")
        j.command(f"mv /tmp/.config ~/")
        j.command(f"gcloud config set account {GCLOUD_USER_ACCOUNT}")
        j.command(f"gcloud config set project {GCLOUD_PROJECT}")

        # run majiq build
        #j.command(f"gsutil -u {GCLOUD_PROJECT} -m cp gs://gtex-resources/GENCODE/gencode.v26.GRCh38.ERCC.genes.collapsed_only.gtf .")
        j.command(f"mv {genes_gtf} gencode.gff3")
        j.command(f"mv {input_read_data.bam} {sample_id}.bam")
        j.command(f"mv {input_read_data.bai} {sample_id}.bam.bai")

        j.command(f"echo '[info]' >> majiq_build.cfg")
        j.command(f"echo 'readlen={metadata_row['read length (rnaseqc)']}' >> majiq_build.cfg")
        j.command(f"echo 'bamdirs=.' >> majiq_build.cfg")
        j.command(f"echo 'genome=hg38' >> majiq_build.cfg")
        j.command(f"echo 'strandness={'None' if metadata_row['stranded? (rnaseqc)'] == 'no' else 'reverse'}' >> majiq_build.cfg")
        j.command(f"echo '[experiments]' >> majiq_build.cfg")
        j.command(f"echo '{sample_id}={sample_id}' >> majiq_build.cfg")

        j.command(f"cat majiq_build.cfg >> {j.logfile}")
        j.command(f"majiq build gencode.gff3 -c majiq_build.cfg -j 1 -o majiq_build_{sample_id} >> {j.logfile}")

        j.command(f"tar czf majiq_build_{sample_id}.tar.gz majiq_build_{sample_id}")
        j.command(f"cp majiq_build_{sample_id}.tar.gz {j.output_tar_gz}")

        #j.command(f"ls -lh . >> {j.logfile}")
        #j.command(f"echo ls majiq_build_{sample_id} >> {j.logfile}")
        #j.command(f"ls -1 majiq_build_{sample_id} >> {j.logfile}")
        j.command(f"echo --- done  {output_file_path} >> {j.logfile}")

        # copy output
        b.write_output(j.output_tar_gz, output_file_path)
        b.write_output(j.logfile, os.path.join(output_dir, f"majiq_build_{sample_id}.log"))

    b.run()

    if isinstance(backend, hb.ServiceBackend):
        backend.close()


if __name__ == "__main__":
    main()
