import argparse
import datetime
import hailtop.batch as hb
import logging
import os

from sample_metadata.utils import get_joined_metadata_df

logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s', level=logging.INFO)
logger = logging.getLogger(__name__)


GCLOUD_PROJECT = "seqr-project"
GCLOUD_USER_ACCOUNT = "weisburd@broadinstitute.org"
GCLOUD_CREDENTIALS_LOCATION = "gs://weisburd-misc/creds"

def main():
    rnaseq_sample_metadata_df = get_joined_metadata_df()

    p = argparse.ArgumentParser()
    grp = p.add_mutually_exclusive_group(required=True)
    grp.add_argument("--cluster", action="store_true", help="Batch: submit to cluster")
    grp.add_argument("--local", action="store_true", help="Batch: run locally")
    p.add_argument("-r", "--raw", action="store_true", help="Batch: run directly on the machine, without using a docker image")
    p.add_argument("--batch-billing-project", default="tgg-rare-disease", help="Batch: billing project. Required if submitting to cluster.")
    p.add_argument("--batch-job-name", help="Batch: (optional) job name")
    p.add_argument("--cpu", type=float, help="Batch: (optional) number of CPUs to request", default=1, choices=[0.25, 0.5, 1, 2, 4, 8, 16])
    p.add_argument("--memory", type=float, help="Batch: (optional) memory in gigabytes", default=3.75)
    p.add_argument("--force", action="store_true", help="Recompute and overwrite cached or previously computed data")
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

    working_dir = "/"

    # see https://hail.zulipchat.com/#narrow/stream/223457-Batch-support/topic/auth.20as.20user.20account for more details
    if args.local:
        if args.raw:
            backend = hb.LocalBackend()
            working_dir = os.path.abspath(os.path.dirname(__name__))
        else:
            backend = hb.LocalBackend(gsa_key_file=os.path.expanduser("~/.config/gcloud/misc-270914-cb9992ec9b25.json"))
    else:
        backend = hb.ServiceBackend(args.batch_billing_project)

    b = hb.Batch(backend=backend, name=args.batch_job_name)

    print("Working dir: " + working_dir)

    # define workflow inputs
    #if args.local:
    #    genes_gtf = b.read_input("gencode.v26.annotation.gff3", extension=".gff3")
    #else:
    #    genes_gtf = b.read_input("gs://macarthurlab-rnaseq/ref/gencode.v26.annotation.GRCh38.gff3", extension=".gff3")

    # define parallel execution for samples
    for sample_id in sample_ids:
        metadata_row = rnaseq_sample_metadata_df.loc[sample_id]
        batch_name = metadata_row['star_pipeline_batch']

        # set job inputs & outputs
        if args.local:
            input_bam, input_bai = [
                os.path.join("/Users/weisburd/project__rnaseq/data/samples/expression/bams/", os.path.basename(path)) for path in (metadata_row['star_bam'], metadata_row['star_bai'])
            ]  # .replace(".bam", ".subset.bam")
            output_dir = "."
        else:
            input_bam, input_bai = metadata_row['star_bam'], metadata_row['star_bai']
            output_dir = f"gs://macarthurlab-rnaseq/{batch_name}/fraser_count_rna/"

            #input_bam = "gs://macarthurlab-rnaseq/temp/MUN_FAM5_SIBLINGMDC1A_01_R1.Aligned.sortedByCoord.out.subset.bam"
            #input_bai = "gs://macarthurlab-rnaseq/temp/MUN_FAM5_SIBLINGMDC1A_01_R1.Aligned.sortedByCoord.out.subset.bam.bai"

        print("Input bam: ", input_bam)
        input_read_data = b.read_input_group(bam=input_bam, bai=input_bai)
        output_file_path = os.path.join(output_dir, f"fraser_count_rna_{sample_id}.tar.gz")

        # check if output file already exists
        if args.local:
            bam_size = None
        else:
            import hail as hl  # used for hadoop file utils
            if hl.hadoop_is_file(output_file_path) and not args.force:
                logger.info(f"{sample_id} output file already exists: {output_file_path}. Skipping...")
                continue

            file_stats = hl.hadoop_stat(metadata_row['star_bam'])
            bam_size = int(round(file_stats['size_bytes']/10.**9))

        # define majiq build commands for this sample
        j = b.new_job(name=args.batch_job_name)
        if not args.raw:
            j.image("weisburd/gagneurlab@sha256:7af235e4beba907569bacbae31271ecb9c6837337b2722da8cfb0587d181d1e8")
        if bam_size:
            j.storage(f'{bam_size*2}Gi')
        j.cpu(args.cpu)  # Batch default is 1
        if args.memory:
            j.memory(f"{args.memory}G")  # Batch default is 3.75G
        else:
            j.memory(f"{3.75*args.cpu}G")  # Batch default is 3.75G
        logger.info(f'Requesting: {j._storage or "default"} storage, {j._cpu or "default"} CPU, {j._memory or "default"} memory')

        # switch to user account
        #j.command(f"gcloud auth activate-service-account --key-file /gsa-key/key.json")
        #j.command(f"gsutil -m cp -r {GCLOUD_CREDENTIALS_LOCATION}/.config /tmp/")
        #j.command(f"mv ~/.config ~/.config."+datetime.datetime.now().strftime("%Y%m%d-%H%M%S"))
        #j.command(f"mv /tmp/.config ~/")
        #j.command(f"gcloud config set account {GCLOUD_USER_ACCOUNT}")
        #j.command(f"gcloud config set project {GCLOUD_PROJECT}")

        #j.command("((while true; do uptime; sleep 30; done) & )")
        j.command("set -x")

        if not args.raw:
            j.command(f"cp {input_read_data.bam} {os.path.join(working_dir, sample_id)}.bam >> {j.logfile}")
            j.command(f"cp {input_read_data.bai} {os.path.join(working_dir, sample_id)}.bam.bai >> {j.logfile}")
            j.command(f"touch {os.path.join(working_dir, sample_id)}.bam.bai >> {j.logfile}")
            bam_path = os.path.join(working_dir, sample_id) + ".bam"
        else:
            bam_path = f"{input_read_data.bam}"

        j.command("bla")
        j.command(f"pwd >> {j.logfile}")
        j.command(f"ls >> {j.logfile}")
        j.command(f"date >> {j.logfile}")
        j.command(f"Rscript --vanilla {os.path.join(working_dir, 'countRNA.R')} --num-threads {args.cpu} {sample_id} {bam_path} >> {j.logfile}")
        j.command(f"ls . >> {j.logfile}")
        j.command(f"tar czf fraser_count_rna_{sample_id}.tar.gz cache")
        j.command(f"cp fraser_count_rna_{sample_id}.tar.gz {j.output_tar_gz}")

        #j.command(f"ls -lh . >> {j.logfile}")
        j.command(f"echo Done: {output_file_path} >> {j.logfile}")
        j.command(f"date >> {j.logfile}")

        # copy output
        b.write_output(j.output_tar_gz, output_file_path)
        print("Output file path: ", output_file_path)

        log_file_path = os.path.join(output_dir, f"fraser_count_rna_{sample_id}.log")
        b.write_output(j.logfile, log_file_path)
        print("Log file path: ", log_file_path)

    b.run()

    if isinstance(backend, hb.ServiceBackend):
        backend.close()


if __name__ == "__main__":
    main()
