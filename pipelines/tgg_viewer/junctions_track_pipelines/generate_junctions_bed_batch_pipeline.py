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

DOCKER_IMAGE = "weisburd/convert-sj-out-tab-to-junctions-bed@sha256:b06670b14608f80679b1e20aea9b35cf5fd19b3f0b12f6839416b366a7cb7ee5"
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
    df.loc[source_df.sample_id, 'star_SJ_out_tab'] = source_df['star_SJ_out_tab']
    df.loc[source_df.sample_id, 'output_dir'] = source_df['star_pipeline_batch'].apply(
        lambda batch_name: f"gs://macarthurlab-rnaseq/{batch_name}/junctions_bed_for_igv_js/")

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

    # Generate samples_df with these columns: sample_id, star_SJ_out_tab, output_dir, batch_name
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
        output_dir = metadata_row['output_dir']

        print("Input file: ", metadata_row['star_SJ_out_tab'])
        output_filename = f"{sample_id}.junctions.bed.gz"
        output_bed_gz_file_path = os.path.join(output_dir, output_filename)

        # check if output file already exists
        if hl.hadoop_is_file(output_bed_gz_file_path) and not args.force:
            logger.info(f"{sample_id} output file already exists: {output_bed_gz_file_path}. Skipping...")
            continue

        j = init_job(b, name=f"tab=>bed: {sample_id}", cpu=args.cpu, memory=args.memory, disk_size=5, switch_to_user_account=True, image=DOCKER_IMAGE)

        j.command(f"gsutil -u {GCLOUD_PROJECT} -m cp {metadata_row['star_SJ_out_tab']} .")
        j.command(f"gsutil -u {GCLOUD_PROJECT} -m cp gs://macarthurlab-rnaseq/ref/gencode.v26.annotation.gff3.gz .")
        j.command(f"pwd && ls && date")

        j.command(f"python3 /convert_SJ_out_tab_to_junctions_bed.py -g gencode.v26.annotation.gff3.gz {os.path.basename(metadata_row['star_SJ_out_tab'])}")
        j.command(f"cp {output_filename} {j.output_bed_gz}")
        j.command(f"cp {output_filename}.tbi {j.output_bed_gz_tbi}")
        j.command(f"echo Done: {output_bed_gz_file_path}")
        j.command(f"date")

        # copy output
        b.write_output(j.output_bed_gz, output_bed_gz_file_path)
        b.write_output(j.output_bed_gz_tbi, f"{output_bed_gz_file_path}.tbi")

        print("Output file path: ", output_bed_gz_file_path)

    b.run()

    if isinstance(backend, hb.ServiceBackend):
        backend.close()


if __name__ == "__main__":
    main()
