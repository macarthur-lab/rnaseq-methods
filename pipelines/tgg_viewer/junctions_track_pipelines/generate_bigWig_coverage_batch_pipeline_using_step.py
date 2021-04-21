import argparse
import logging
import pandas as pd
import os

from batch import batch_wrapper

logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s', level=logging.INFO)
logger = logging.getLogger(__name__)

DOCKER_IMAGE = "gcr.io/broad-cga-francois-gtex/gtex_rnaseq:V9"
GCLOUD_PROJECT = "seqr-project"
GCLOUD_USER_ACCOUNT = "weisburd@broadinstitute.org"
GCLOUD_CREDENTIALS_LOCATION = "gs://weisburd-misc/creds"


def get_sample_table(parser: argparse.ArgumentParser, expected_columns: tuple) -> pd.DataFrame:
    parser.add_argument("sample_table", help=f"Input .tsv with one row per sample, and columns {', '.join(expected_columns)}")
    args = parser.parse_known_args()  # get the sample table

    df = pd.read_table(args.sample_table)
    missing_columns = set(expected_columns) - set(df.columns)
    if missing_columns:
        parser.error(f"{args.sample_table} needs to have these columns: " + ", ".join(missing_columns))

    logging.info(f"Parsed {len(df)} rows from {args.sample_table}")

    grp = parser.add_mutually_exclusive_group()
    grp.add_argument("-n", "--num-samples", action="append", help="num samples to process")
    grp.add_argument("-s", "--sample-id", action="append", help="sample IDs to process (eg. -s sample_id1 -s sample_id2)", choices=set(df['sample_id']))
    args = parser.parse_known_args()

    if args.sample_id or args.num_samples:
        if args.sample_id:
            filtered_df = df[df["sample_id"].isin(args.sample_id)]
        elif args.num_samples:
            filtered_df = df.iloc[:args.num_samples]

        logger.info(f"Keeping {len(filtered_df)} out of {len(df)} rows")
        df = filtered_df

    return df


def main():
    #p = batch_utils.init_arg_parser(default_cpu=4, gsa_key_file=os.path.expanduser("~/.config/gcloud/misc-270914-cb9992ec9b25.json"))

    with batch_wrapper.batch_pipeline() as bp:
        p = bp.get_argument_parser()

        df = get_sample_table(p, expected_columns=("sample_id", "bam_path", "bai_path", "output_dir"))

        # set job inputs & outputs

        s = bp.new_step()

        for sample_id in df.sample_id:
            metadata_row = df.loc[sample_id]

            output_filename = f"{sample_id}.bigWig"
            output_file_path = os.path.join(output_dir, output_filename)

            s.input(input_bam)
            s.input(input_bai)
            s.input("gs://gtex-resources/references/GRCh38.chrsizes")
            s.command(f"python3 /src/bam2coverage.py {sample_id}.bam GRCh38.chrsizes {sample_id}")
            s.output(output_file_path)

            # check if output file already exists
            #if hl.hadoop_is_file(output_file_path) and not args.force:
            #    logger.info(f"{sample_id} output file already exists: {output_file_path}. Skipping...")
            #    continue

            #file_stats = hl.hadoop_stat(metadata_row['bam_path'])
            #bam_size = int(round(file_stats['size_bytes']/10.**9))
            #disk_size = bam_size * 2

            #j = batch_utils.init_job(batch, f"bam=>bigWig: {sample_id}", cpu=args.cpu, memory=args.memory, disk_size=disk_size, image=DOCKER_IMAGE)
            #batch_utils.switch_gcloud_auth_to_user_account(j, GCLOUD_CREDENTIALS_LOCATION, GCLOUD_USER_ACCOUNT, GCLOUD_PROJECT)

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
            batch.write_output(j.output_bigWig, output_file_path)

            print("Output file path: ", output_file_path)


if __name__ == "__main__":
    main()
