import hail as hl  # used for hadoop file utils
import logging
import os
import pandas as pd

from batch import batch_utils
from sample_metadata.rnaseq_metadata_utils import get_rnaseq_metadata_joined_with_paths_df

hl.init(log="/dev/null")

logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s', level=logging.INFO)
logger = logging.getLogger(__name__)

DOCKER_IMAGE = "weisburd/portcullis@sha256:430c24f61f8c580854857f655e7f460494414403b5d15a2963840cab0ba334c2"
GCLOUD_PROJECT = "seqr-project"
GCLOUD_USER_ACCOUNT = "weisburd@broadinstitute.org"
GCLOUD_CREDENTIALS_LOCATION = "gs://weisburd-misc/creds"


def main():
    metadata_df = get_rnaseq_metadata_joined_with_paths_df()
    tissues = set([b.strip() for b in metadata_df["imputed tissue"] if isinstance(b, str)])
    star_pipeline_batches = set([b for b in metadata_df["star_pipeline_batch"] if b])

    p = batch_utils.init_arg_parser(gsa_key_file=os.path.expanduser("~/.config/gcloud/misc-270914-cb9992ec9b25.json"))
    p.add_argument("tissue_or_sample_id", nargs="+", choices={"all",} | tissues | star_pipeline_batches | set(metadata_df['sample_id']))
    p.add_argument("-gtex", "--include-gtex-samples", action="store_true")
    args = p.parse_args()

    if not args.include_gtex_samples:
        metadata_df = metadata_df[~metadata_df.sample_id.str.lower().str.startswith("gtex")]

    # Generate samples_df with these columns: sample_id, bam_path, bai_path, output_dir, batch_name, sex, RIN, ancestry, etc.
    all_samples_df = None
    for name in args.tissue_or_sample_id:
        if name == "all":
            samples_df = metadata_df
        elif name in star_pipeline_batches:
            samples_df = metadata_df[metadata_df['star_pipeline_batch'] == name]
        elif name in tissues:
            samples_df = metadata_df[metadata_df['imputed tissue'] == name]
        elif name in set(metadata_df.sample_id):
            samples_df = metadata_df[metadata_df.sample_id == name]
        else:
            p.error(f"Unexpected name {name}")

    if all_samples_df is None:
        all_samples_df = samples_df
    else:
        all_samples_df = pd.concat([all_samples_df, samples_df], axis="rows")

    all_samples_df.loc[:, 'output_dir'] = all_samples_df['star_pipeline_batch'].apply(lambda batch_name: f"gs://tgg-rnaseq/{batch_name}/portcullis/")
    all_samples_df = all_samples_df.set_index('sample_id', drop=False)
    samples_df = all_samples_df
    logger.info(f"Processing {len(samples_df)} sample ids: {', '.join(samples_df.sample_id[:20])}")

    # see https://hail.zulipchat.com/#narrow/stream/223457-Batch-support/topic/auth.20as.20user.20account for more details
    with batch_utils.run_batch(args) as batch:
        for sample_id in samples_df.sample_id:
            metadata_row = samples_df.loc[sample_id]

            # set job inputs & outputs
            input_bam, input_bai = metadata_row['star_bam'], metadata_row['star_bai']
            output_dir = metadata_row['output_dir']

            print("Input bam: ", input_bam)
            output_all_junctions_tab_filename = f"{sample_id}.portcullis_all.junctions.tab.gz"
            output_all_junctions_tab_file_path = os.path.join(output_dir, output_all_junctions_tab_filename)
            output_filtered_junctions_tab_filename = f"{sample_id}.portcullis_filtered.pass.junctions.tab.gz"
            output_filtered_junctions_tab_file_path = os.path.join(output_dir, output_filtered_junctions_tab_filename)
            output_bed_filename = f"{sample_id}.portcullis_filtered.pass.junctions.bed.gz"
            output_bed_file_path = os.path.join(output_dir, output_bed_filename)

            # check if output file already exists
            if hl.hadoop_is_file(output_filtered_junctions_tab_file_path) and hl.hadoop_is_file(output_bed_file_path) and not args.force:
                logger.info(f"{sample_id} output files already exist: {output_filtered_junctions_tab_file_path} {output_bed_file_path}. Skipping...")
                continue

            file_stats = hl.hadoop_stat(metadata_row['star_bam'])
            bam_size = int(round(file_stats['size_bytes']/10.**9))
            disk_size = bam_size * 1.5

            j = batch_utils.init_job(batch, f"portcullis: {sample_id}", cpu=8, disk_size=disk_size, image=DOCKER_IMAGE)
            batch_utils.switch_gcloud_auth_to_user_account(j, GCLOUD_CREDENTIALS_LOCATION, GCLOUD_USER_ACCOUNT, GCLOUD_PROJECT)

            local_fasta = batch_utils.localize_file(j, batch_utils.HG38_REF_PATHS.fasta, use_gcsfuse=True)
            local_bam = batch_utils.localize_file(j, input_bam, use_gcsfuse=True)
            local_bai = batch_utils.localize_file(j, input_bai, use_gcsfuse=True)
            local_gencode_gtf = batch_utils.localize_file(j, batch_utils.HG38_REF_PATHS.gencode_v36_gtf, use_gcsfuse=True)

            #j.command(f"touch {local_bai}")

            j.command(f"pwd && ls && date")

            if metadata_row["stranded? (rnaseqc)"] == "yes":
                stranded_arg = "--strandedness unstranded"
            elif metadata_row["stranded? (rnaseqc)"] == "no":
                stranded_arg = "--strandedness firststrand"
            else:
                raise ValueError(f"Unexpected 'stranded' value for {metadata_row['sample_id']} : {metadata_row['stranded? (rnaseqc)']}")

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
