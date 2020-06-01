import argparse
import datetime
import hailtop.batch as hb
import hashlib
import logging
import os

from sample_metadata.utils import get_joined_metadata_df

logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s', level=logging.INFO)
logger = logging.getLogger(__name__)


GCLOUD_PROJECT = "seqr-project"
GCLOUD_USER_ACCOUNT = "weisburd@broadinstitute.org"
GCLOUD_CREDENTIALS_LOCATION = "gs://weisburd-misc/creds"


def init_job(
    batch,
    name: str = None,
    image: str = None,
    disk_size: float = None,
    cpu: float = 1,
    memory: float = None,
    switch_to_user_account: bool = False,
):

    j = batch.new_job(name=name)
    if image:
        j.image(image)

    if disk_size:
        j.storage(f'{disk_size}Gi')

    j.cpu(cpu)  # Batch default is 1
    if memory:
        j.memory(f"{memory}G")  # Batch default is 3.75G
    else:
        j.memory(f"{3.75*cpu}G")  # Batch default is 3.75G

    logger.info(f'Requesting: {j._storage or "default"} storage, {j._cpu or "default"} CPU, {j._memory or "default"} memory')

    # switch to user account
    if switch_to_user_account:
        j.command(f"gcloud auth activate-service-account --key-file /gsa-key/key.json")
        j.command(f"gsutil -m cp -r {GCLOUD_CREDENTIALS_LOCATION}/.config /tmp/")
        j.command(f"mv ~/.config ~/.config."+datetime.datetime.now().strftime("%Y%m%d-%H%M%S"))
        j.command(f"mv /tmp/.config ~/")
        j.command(f"gcloud config set account {GCLOUD_USER_ACCOUNT}")
        j.command(f"gcloud config set project {GCLOUD_PROJECT}")

    #j.command("((while true; do uptime; sleep 30; done) & )")
    j.command("set -x")

    return j


def get_batch_label(sample_ids):
    byte_string = ", ".join(sorted(sample_ids)).encode()
    h = hashlib.md5(byte_string).hexdigest().upper()
    return f"{len(sample_ids)}_samples_{h[:10]}"


def main():
    rnaseq_sample_metadata_df = get_joined_metadata_df()

    p = argparse.ArgumentParser()
    grp = p.add_mutually_exclusive_group(required=True)
    grp.add_argument("--cluster", action="store_true", help="Batch: submit to cluster")
    grp.add_argument("--local", action="store_true", help="Batch: run locally")
    p.add_argument("-r", "--raw", action="store_true", help="Batch: run directly on the machine, without using a docker image")
    p.add_argument("--batch-billing-project", default="tgg-rare-disease", help="Batch: billing project. Required if submitting to cluster.")
    p.add_argument("--batch-job-name", help="Batch: (optional) job name")
    p.add_argument("--batch-temp-bucket", default="macarthurlab-rnaseq", help="Batch: bucket where it stores temp files. "
        "The batch service-account must have Admin permissions for this bucket. These can be added by running "
        "gsutil iam ch serviceAccount:[SERVICE_ACCOUNT_NAME]:objectAdmin gs://[BUCKET_NAME]")
    p.add_argument("-t", "--cpu", type=float, help="Batch: (optional) number of CPUs (eg. 0.5)", default=1, choices=[0.25, 0.5, 1, 2, 4, 8, 16])
    p.add_argument("-m", "--memory", type=float, help="Batch: (optional) memory in gigabytes (eg. 3.75)", default=3.75)
    p.add_argument("--force", action="store_true", help="Recompute and overwrite cached or previously computed data")
    grp = p.add_mutually_exclusive_group(required=True)
    grp.add_argument("-b", "--rnaseq-batch-name", nargs="*", help="RNA-seq batch names to process (eg. -b batch1 batch2)", choices=set(rnaseq_sample_metadata_df['star_pipeline_batch']))
    grp.add_argument("-s", "--rnaseq-sample-id", nargs="*", help="RNA-seq sample IDs to process (eg. -s sample1 sample2)", choices=set(rnaseq_sample_metadata_df['sample_id']))
    args = p.parse_args()

    #logger.info("\n".join(df.columns))

    if args.rnaseq_batch_name:
        batch_names = args.rnaseq_batch_name
        sample_ids = rnaseq_sample_metadata_df[rnaseq_sample_metadata_df['star_pipeline_batch'].isin(batch_names)].sample_id
    elif args.rnaseq_sample_id:
        sample_ids = args.rnaseq_sample_id
    else:
        p.error("Must specify -b or -s")

    logger.info(f"Processing {len(sample_ids)} sample ids: {', '.join(sample_ids)}")

    working_dir = "/"

    # see https://hail.zulipchat.com/#narrow/stream/223457-Batch-support/topic/auth.20as.20user.20account for more details
    if args.local:
        if args.raw:
            backend = hb.LocalBackend()
            working_dir = "./"
        else:
            backend = hb.LocalBackend(gsa_key_file=os.path.expanduser("~/.config/gcloud/misc-270914-cb9992ec9b25.json"))
    else:
        backend = hb.ServiceBackend(billing_project=args.batch_billing_project) #, bucket=args.batch_temp_bucket)

    b = hb.Batch(backend=backend, name=args.batch_job_name)

    print("Working dir: " + working_dir)

    # define workflow inputs
    #if args.local:
    #    genes_gtf = b.read_input("gencode.v26.annotation.gff3", extension=".gff3")
    #else:
    #    genes_gtf = b.read_input("gs://macarthurlab-rnaseq/ref/gencode.v26.annotation.GRCh38.gff3", extension=".gff3")

    split_reads_samples = []
    split_reads_output_files = []
    split_reads_jobs = []
    j_extract_splice_junctions = init_job(
        b, name=f"Extract Splice-junctions {args.batch_job_name or ''}", switch_to_user_account=False,
        image=None if args.raw else "weisburd/gagneurlab@sha256:8e5038eed349438b0c0778ada7bc41b1852e1be3cd69612ae26b28a43619b21d")

    for step in 1, 2:
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
            if step == 1:
                output_file_path = os.path.join(output_dir, f"fraser_count_split_reads_{sample_id}.tar.gz")
            else:
                output_file_path = os.path.join(output_dir, f"fraser_count_non_split_reads_{sample_id}.tar.gz")

            # check if output file already exists
            if args.local:
                disk_size = None
            else:
                import hail as hl  # used for hadoop file utils
                if hl.hadoop_is_file(output_file_path) and not args.force:
                    logger.info(f"{sample_id} output file already exists: {output_file_path}. Skipping...")
                    continue

                file_stats = hl.hadoop_stat(metadata_row['star_bam'])
                bam_size = int(round(file_stats['size_bytes']/10.**9))
                disk_size = bam_size * 2

            job_label = f"Count {'Split' if step == 1 else 'Non-split'} Reads"
            j = init_job(b, name=f"{job_label} {args.batch_job_name or ''}", cpu=args.cpu, memory=args.memory, disk_size=disk_size, switch_to_user_account=False,
                                image=None if args.raw else "weisburd/gagneurlab@sha256:8e5038eed349438b0c0778ada7bc41b1852e1be3cd69612ae26b28a43619b21d")

            if not args.raw:
                j.command(f"cp {input_read_data.bam} {os.path.join(working_dir, sample_id)}.bam")
                j.command(f"cp {input_read_data.bai} {os.path.join(working_dir, sample_id)}.bam.bai")
                j.command(f"touch {os.path.join(working_dir, sample_id)}.bam.bai")
                bam_path = os.path.join(working_dir, sample_id) + ".bam"
            else:
                bam_path = f"{input_read_data.bam}"

            j.command(f"pwd && ls && date")
            if step == 1:
                script = os.path.join(working_dir, 'countSplitReads.R')
                j.command(f"Rscript --vanilla {script} {sample_id} {bam_path}")
            else:
                script = os.path.join(working_dir, 'countNonSplitReads.R')
                j.command(f"Rscript --vanilla {script} {sample_id} {bam_path} {j_extract_splice_junctions.splice_junctions_RDS}")

            j.command(f"ls .")
            j.command(f"tar czf {j.output_tar_gz} cache")

            #j.command(f"ls -lh .")
            j.command(f"echo Done: {output_file_path}")
            j.command(f"date")

            # copy output
            b.write_output(j.output_tar_gz, output_file_path)
            print("Output file path: ", output_file_path)

            if step == 1:
                split_reads_samples.append(sample_id)
                split_reads_output_files.append(output_file_path)
                split_reads_jobs.append(j)

        if step == 1:
            batch_label = get_batch_label(split_reads_samples)
            for j in split_reads_jobs:
                j_extract_splice_junctions.depends_on(j)

            j_extract_splice_junctions.command(f"gsutil -m cp {' '.join(split_reads_output_files)} .")
            j_extract_splice_junctions.command(f"for i in fraser_count_split_reads*.tar.gz; do tar xzf $i; done")
            j_extract_splice_junctions.command(f"pwd && ls && date")
            j_extract_splice_junctions.command(f"Rscript --vanilla {os.path.join(working_dir, 'extractSpliceJunctions.R')}")
            j_extract_splice_junctions.command(f"ls .")
            j_extract_splice_junctions.command(f"cp splice_junctions.RDS {j_extract_splice_junctions.splice_junctions_RDS}")
            output_file_path = os.path.join("gs://macarthurlab-rnaseq", f"splice_junctions_{batch_label}.RDS")

            b.write_output(j_extract_splice_junctions.splice_junctions_RDS, output_file_path)
            print("Output file path: ", output_file_path)

    b.run()

    if isinstance(backend, hb.ServiceBackend):
        backend.close()


if __name__ == "__main__":
    main()
