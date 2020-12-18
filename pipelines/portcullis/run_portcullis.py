import hail as hl  # used for hadoop file utils
import logging
import os
import pandas as pd

from batch import batch_utils
from sample_metadata.rnaseq_metadata_utils import get_joined_metadata_df, get_gtex_rnaseq_sample_metadata_df

hl.init(log="/dev/null")


logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s', level=logging.INFO)
logger = logging.getLogger(__name__)

DOCKER_IMAGE = "weisburd/portcullis@sha256:430c24f61f8c580854857f655e7f460494414403b5d15a2963840cab0ba334c2"
GCLOUD_PROJECT = "seqr-project"
GCLOUD_USER_ACCOUNT = "weisburd@broadinstitute.org"
GCLOUD_CREDENTIALS_LOCATION = "gs://weisburd-misc/creds"


def transfer_metadata_columns_from_df(samples_df, source_df):
    df = pd.DataFrame()
    df.loc[:, 'sample_id'] = source_df['sample_id']
    df.loc[source_df.sample_id, 'stranded'] = source_df['stranded? (rnaseqc)']
    df.loc[source_df.sample_id, 'batch'] = source_df['star_pipeline_batch']
    df.loc[source_df.sample_id, 'bam_path'] = source_df['star_bam']
    df.loc[source_df.sample_id, 'bai_path'] = source_df['star_bai']

    df.loc[source_df.sample_id, 'output_dir'] = source_df['star_pipeline_batch'].apply(
        lambda batch_name: f"gs://macarthurlab-rnaseq/{batch_name}/portcullis/")

    return pd.concat([samples_df, df], axis="rows")


def main():
    rnaseq_sample_metadata_df = get_joined_metadata_df()
    #gtex_rnaseq_sample_metadata_df = get_gtex_rnaseq_sample_metadata_df()

    p = batch_utils.init_arg_parser(gsa_key_file=os.path.expanduser("~/.config/gcloud/misc-270914-cb9992ec9b25.json"))
    grp = p.add_mutually_exclusive_group(required=True)
    grp.add_argument("-b", "--rnaseq-batch-name", nargs="*", help="RNA-seq batch names to process (eg. -b batch1 batch2)",
                     choices=set(rnaseq_sample_metadata_df['star_pipeline_batch']))
    grp.add_argument("-s", "--rnaseq-sample-id", nargs="*", help="RNA-seq sample IDs to process (eg. -s sample1 sample2)",
                     choices=set(rnaseq_sample_metadata_df['sample_id']))
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
    with batch_utils.run_batch(args) as batch:
        for sample_id in samples_df.sample_id:
            metadata_row = samples_df.loc[sample_id]

            # set job inputs & outputs
            input_bam, input_bai = metadata_row['bam_path'], metadata_row['bai_path']
            output_dir = metadata_row['output_dir']

            print("Input bam: ", input_bam)
            output_all_junctions_tab_filename = f"{sample_id}.portcullis_all.junctions.tab.gz"
            output_all_junctions_tab_file_path = os.path.join(output_dir, output_all_junctions_tab_filename)
            output_filtered_junctions_tab_filename = f"{sample_id}.portcullis_filtered.pass.junctions.tab.gz"
            output_filtered_junctions_tab_file_path = os.path.join(output_dir, output_filtered_junctions_tab_filename)
            output_bed_filename = f"{sample_id}.portcullis_filtered.pass.junctions.bed.gz"
            output_bed_file_path = os.path.join(output_dir, output_bed_filename)

            # check if output file already exists
            if hl.hadoop_is_file(output_filtered_junctions_tab_file_path) and hl.hadoop_isfile(output_bed_file_path) and not args.force:
                logger.info(f"{sample_id} output files already exist: {output_filtered_junctions_tab_file_path} {output_bed_file_path}. Skipping...")
                continue

            file_stats = hl.hadoop_stat(metadata_row['bam_path'])
            bam_size = int(round(file_stats['size_bytes']/10.**9))
            disk_size = bam_size * 1.5

            j = batch_utils.init_job(batch, f"portcullis: {sample_id}", cpu=4, disk_size=disk_size, image=DOCKER_IMAGE)
            batch_utils.switch_gcloud_auth_to_user_account(j, GCLOUD_CREDENTIALS_LOCATION, GCLOUD_USER_ACCOUNT, GCLOUD_PROJECT)

            local_fasta = batch_utils.localize_file(j, batch_utils.HG38_REF_PATHS.fasta, use_gcsfuse=True)
            local_bam = batch_utils.localize_file(j, input_bam, use_gcsfuse=True)
            local_bai = batch_utils.localize_file(j, input_bai, use_gcsfuse=True)
            local_gencode_gtf = batch_utils.localize_file(j, batch_utils.HG38_REF_PATHS.gencode_v36_gtf, use_gcsfuse=True)

            #j.command(f"touch {local_bai}")

            j.command(f"pwd && ls && date")

            if metadata_row["stranded"] == "yes":
                stranded_arg = "--strandedness unstranded"
            elif metadata_row["stranded"] == "no":
                stranded_arg = "--strandedness firststrand"
            else:
                raise ValueError(f"Unexpected 'stranded' value for {metadata_row['sample_id']} : {metadata_row['stranded']}")

            j.command(f"time portcullis full -t 4 --orientation FR {stranded_arg} -r {local_gencode_gtf} -v {local_fasta} {local_bam}")
            #j.command(f"find .")
            
            j.command(f"gzip -c ./portcullis_out/2-junc/portcullis_all.junctions.tab > {output_all_junctions_tab_filename}")
            j.command(f"gzip -c ./portcullis_out/3-filt/portcullis_filtered.pass.junctions.tab > {output_filtered_junctions_tab_filename}")
            j.command(f"gzip -c ./portcullis_out/3-filt/portcullis_filtered.pass.junctions.bed > {output_bed_filename}")
            #j.command(f"tabix {output_bed_filename}")  # need to sort the bed file first
            j.command(f"gsutil -m cp {output_all_junctions_tab_filename} {output_all_junctions_tab_file_path}")
            j.command(f"gsutil -m cp {output_filtered_junctions_tab_filename} {output_filtered_junctions_tab_file_path}")
            j.command(f"gsutil -m cp {output_bed_filename} {output_bed_file_path}")
            #j.command(f"gsutil -m cp {output_bed_filename}.tbi {output_bed_file_path}.tbi")
            j.command(f"echo Done: {output_filtered_junctions_tab_file_path}")
            j.command(f"date")

            print("Output file path: ", output_filtered_junctions_tab_file_path)


if __name__ == "__main__":
    main()
