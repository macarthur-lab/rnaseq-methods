import argparse
import datetime
import hail as hl  # used for hadoop file utils
import hailtop.batch as hb
import hashlib
import logging
import os
import pandas as pd
import sys

from sample_metadata.utils import get_joined_metadata_df, get_gtex_rnaseq_sample_metadata_df

logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s', level=logging.INFO)
logger = logging.getLogger(__name__)

DOCKER_IMAGE = "gcr.io/broad-cga-francois-gtex/gtex_rnaseq:V9"
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

    j.command("set -x")

    return j


def transfer_metadata_columns_from_df(samples_df, source_df):
    df = pd.DataFrame()
    df.loc[:, 'sample_id'] = source_df['sample_id']
    df.loc[source_df.sample_id, 'batch'] = source_df['star_pipeline_batch']
    df.loc[source_df.sample_id, 'bam_path'] = source_df['star_bam']
    df.loc[source_df.sample_id, 'bai_path'] = source_df['star_bai']

    df.loc[source_df.sample_id, 'output_dir'] = source_df['star_pipeline_batch'].apply(
        lambda batch_name: f"gs://macarthurlab-rnaseq/{batch_name}/bigWig/")

    return pd.concat([samples_df, df], axis="rows")


def main():
    rnaseq_sample_metadata_df = get_joined_metadata_df()
    #gtex_rnaseq_sample_metadata_df = get_gtex_rnaseq_sample_metadata_df()

    p = argparse.ArgumentParser()
    p.add_argument("--batch-billing-project", default="tgg-rare-disease", help="Batch: billing project. Required if submitting to cluster.")
    p.add_argument("--batch-name", help="Batch: (optional) batch name")
    p.add_argument("--batch-temp-bucket", default="macarthurlab-rnaseq", help="Batch: bucket where it stores temp files. "
        "The batch service-account must have Admin permissions for this bucket. These can be added by running "
        "gsutil iam ch serviceAccount:[SERVICE_ACCOUNT_NAME]:objectAdmin gs://[BUCKET_NAME]")
    p.add_argument("-t", "--cpu", type=float, help="Batch: (optional) number of CPUs (eg. 0.5)", default=1, choices=[0.25, 0.5, 1, 2, 4, 8, 16])
    p.add_argument("-m", "--memory", type=float, help="Batch: (optional) memory in gigabytes (eg. 3.75)", default=3.75)
    p.add_argument("--force", action="store_true", help="Recompute and overwrite cached or previously computed data")
    grp = p.add_mutually_exclusive_group(required=True)
    grp.add_argument("-b", "--rnaseq-batch-name", nargs="*", help="RNA-seq batch names to process (eg. -b batch1 batch2)",
        choices=set(rnaseq_sample_metadata_df['star_pipeline_batch']) | set(["gtex_muscle", "gtex_fibroblasts", "gtex_blood"]))
    grp.add_argument("-s", "--rnaseq-sample-id", nargs="*", help="RNA-seq sample IDs to process (eg. -s sample1 sample2)",
        choices=set(rnaseq_sample_metadata_df['sample_id']) | set(['GTEX-1LG7Z-0005-SM-DKPQ6', 'GTEX-PX3G-0006-SM-5SI7E', 'GTEX-1KXAM-0005-SM-DIPEC']))
    args = p.parse_args()

    # Generate samples_df with these columns: sample_id, bam_path, bai_path, output_dir, batch_name, sex, RIN, ancestry, etc.
    samples_df = pd.DataFrame()
    if args.rnaseq_batch_name:
        for batch_name in args.rnaseq_batch_name:
            df = rnaseq_sample_metadata_df[rnaseq_sample_metadata_df['star_pipeline_batch'] == batch_name]
            samples_df = transfer_metadata_columns_from_df(samples_df, df)

    elif args.rnaseq_sample_id:
        df = rnaseq_sample_metadata_df[rnaseq_sample_metadata_df.sample_id.isin(set(args.rnaseq_sample_id))]
        samples_df = transfer_metadata_columns_from_df(samples_df, df)
    else:
        p.error("Must specify -b or -s")

    logger.info(f"Processing {len(samples_df)} sample ids: {', '.join(samples_df.sample_id[:20])}")

    # see https://hail.zulipchat.com/#narrow/stream/223457-Batch-support/topic/auth.20as.20user.20account for more details
    backend = hb.ServiceBackend(billing_project=args.batch_billing_project, bucket=args.batch_temp_bucket)
    b = hb.Batch(backend=backend, name=args.batch_name)
    for sample_id in samples_df.sample_id:
        metadata_row = samples_df.loc[sample_id]

        # set job inputs & outputs
        input_bam, input_bai = metadata_row['bam_path'], metadata_row['bai_path']
        output_dir = metadata_row['output_dir']

        print("Input bam: ", input_bam)
        output_filename = f"{sample_id}.bigWig"
        output_file_path = os.path.join(output_dir, output_filename)

        # check if output file already exists
        if hl.hadoop_is_file(output_file_path) and not args.force:
            logger.info(f"{sample_id} output file already exists: {output_file_path}. Skipping...")
            continue

        file_stats = hl.hadoop_stat(metadata_row['bam_path'])
        bam_size = int(round(file_stats['size_bytes']/10.**9))
        disk_size = bam_size * 2

        j = init_job(b, name=f"bam=>bigWig: {sample_id}", cpu=args.cpu, memory=args.memory, disk_size=disk_size, switch_to_user_account=True, image=DOCKER_IMAGE)

        j.command(f"gsutil -u {GCLOUD_PROJECT} -m cp {input_bam} {sample_id}.bam")
        j.command(f"gsutil -u {GCLOUD_PROJECT} -m cp {input_bai} {sample_id}.bam.bai")
        j.command(f"gsutil -u {GCLOUD_PROJECT} -m cp gs://gtex-resources/references/GRCh38.chrsizes .")
        j.command(f"touch {sample_id}.bam.bai")

        j.command(f"pwd && ls && date")

        j.command(f"python3 /src/bam2coverage.py {sample_id}.bam GRCh38.chrsizes {sample_id}")
        j.command(f"cp {output_filename} {j.output_bigWig}")
        j.command(f"echo Done: {output_file_path}")
        j.command(f"date")

        # copy output
        b.write_output(j.output_bigWig, output_file_path)

        print("Output file path: ", output_file_path)

    b.run()

    if isinstance(backend, hb.ServiceBackend):
        backend.close()


if __name__ == "__main__":
    main()
