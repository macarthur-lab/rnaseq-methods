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

from batch import batch_utils


logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s', level=logging.INFO)
logger = logging.getLogger(__name__)

DOCKER_IMAGE = "weisburd/convert-sj-out-tab-to-junctions-bed@sha256:33b0ff3334a6f73613c2cf4f44e4a53e95793590c9fd0a92f03c74702feb7595"
GCLOUD_PROJECT = "seqr-project"
GCLOUD_USER_ACCOUNT = "weisburd@broadinstitute.org"
GCLOUD_CREDENTIALS_LOCATION = "gs://weisburd-misc/creds"


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

    p = batch_utils.init_arg_parser(default_cpu=4, gsa_key_file=os.path.expanduser("~/.config/gcloud/misc-270914-cb9992ec9b25.json"))
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
    with batch_utils.run_batch(args) as batch:
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

            j = batch_utils.init_job(batch, name=f"tab=>bed: {sample_id}", cpu=args.cpu, memory=args.memory, disk_size=5, image=DOCKER_IMAGE)
            batch_utils.switch_gcloud_auth_to_user_account(j, GCLOUD_CREDENTIALS_LOCATION, GCLOUD_USER_ACCOUNT, GCLOUD_PROJECT)

            j.command(f"gsutil -u {GCLOUD_PROJECT} -m cp {metadata_row['star_SJ_out_tab']} .")
            j.command(f"gsutil -u {GCLOUD_PROJECT} -m cp gs://macarthurlab-rnaseq/ref/gencode.v26.annotation.gff3.gz .")
            j.command(f"pwd && ls && date")

            j.command(f"python3 /convert_SJ_out_tab_to_junctions_bed.py -g gencode.v26.annotation.gff3.gz {os.path.basename(metadata_row['star_SJ_out_tab'])}")
            j.command(f"cp {output_filename} {j.output_bed_gz}")
            j.command(f"cp {output_filename}.tbi {j.output_bed_gz_tbi}")
            j.command(f"echo Done: {output_bed_gz_file_path}")
            j.command(f"date")

            # copy output
            batch.write_output(j.output_bed_gz, output_bed_gz_file_path)
            batch.write_output(j.output_bed_gz_tbi, f"{output_bed_gz_file_path}.tbi")

            print("Output file path: ", output_bed_gz_file_path)


if __name__ == "__main__":
    main()
