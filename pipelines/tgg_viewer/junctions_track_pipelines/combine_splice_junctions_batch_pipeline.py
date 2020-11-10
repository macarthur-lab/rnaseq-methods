import hail as hl
import logging
import os
from pprint import pprint

from batch import batch_utils
from sample_metadata.rnaseq_metadata_utils import get_joined_metadata_df, get_gtex_rnaseq_sample_metadata_df, \
    ANALYSIS_BATCHES

logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s', level=logging.INFO)
logger = logging.getLogger(__name__)


GCLOUD_PROJECT = "seqr-project"
GCLOUD_USER_ACCOUNT = "weisburd@broadinstitute.org"
GCLOUD_CREDENTIALS_LOCATION = "gs://weisburd-misc/creds"

DOCKER_IMAGE = "weisburd/convert-sj-out-tab-to-junctions-bed@sha256:3be35e1eb34cbea3a1ca32d8108b969793f41071730bc6135366add601c0cc7c"


def combine_splice_junctions(args, batch, batch_name, SJ_out_tab_paths, save_individual_tables, normalize_read_counts, output_dir):

    output_filename = f"combined.{batch_name}.{len(SJ_out_tab_paths)}_samples.SJ.out.tsv"
    output_path = os.path.join(output_dir, output_filename)
    output_path_exists = hl.hadoop_is_file(output_path)

    output_filename2 = output_filename.replace(".SJ.out.tsv", ".junctions.bed.gz")
    output_path2 = os.path.join(output_dir, output_filename2)
    output_path2_exists = hl.hadoop_is_file(output_path2)

    if not args.force and output_path_exists and output_path2_exists:
        logger.info(f"Output files \n{output_path} and \n{output_path2} exist.\n Skipping...")
        return

    j = batch_utils.init_job(batch, f"combine junctions: {batch_name} ({len(SJ_out_tab_paths)} files)", DOCKER_IMAGE if not args.raw else None, args.cpu, args.cpu*3.75)
    batch_utils.switch_gcloud_auth_to_user_account(j, GCLOUD_CREDENTIALS_LOCATION, GCLOUD_USER_ACCOUNT, GCLOUD_PROJECT)

    if args.force or not output_path_exists:
        local_SJ_out_tab_paths = []
        for SJ_out_tab_path in SJ_out_tab_paths:
            local_path = batch_utils.localize_file(j, SJ_out_tab_path, use_gcsfuse=False)
            local_SJ_out_tab_paths.append(local_path)

        #local_gencode_gff_path = batch_utils.localize_file(j, "gencode.v26.annotation.gff3.gz", use_gcsfuse=False)

        save_individual_tables_option = "--save-individual-tables" if save_individual_tables else ""
        normalize_read_counts_option = "--normalize-read-counts" if normalize_read_counts else ""

        local_SJ_out_tab_paths = " ".join(local_SJ_out_tab_paths)
        j.command(f"python3 -u combine_splice_junctions_using_pandas.py -n 20 "
            f"{save_individual_tables_option} "
            f"{normalize_read_counts_option} "
            f"{local_SJ_out_tab_paths}")
        j.command(f"mv combined.{len(SJ_out_tab_paths)}_samples.SJ.out.tsv {output_filename}")
        j.command(f"""gsutil -m cp {output_filename} {output_dir}""")
        input_path_for_step2 = output_filename
    else:
        logger.info(f"Output file {output_path} exists. Skipping...")
        local_path = batch_utils.localize_file(j, output_path, use_gcsfuse=False)
        j.command(f"mv {local_path} .")
        input_path_for_step2 = os.path.basename(local_path)

    if args.force or not output_path2_exists:
        j.command(f"python3 -u convert_SJ_out_tab_to_junctions_bed.py "
            f"-g gencode.v35.annotation.gff3.gz "
            f"{input_path_for_step2}")
        j.command(f"""gsutil -m cp {output_filename2}* {output_dir}""")
    else:
        logger.info(f"Output file {output_path2} exists. Skipping...")

    logger.info(f"Output: {output_path}")
    logger.info(f"Output2: {output_path2}")


if __name__ == "__main__":
    # columns: sample_id, star_pipeline_batch, star_SJ_out_tab, 'imputed sex', 'imputed tissue', 'stranded? (rnaseqc)', 'read length (rnaseqc)'
    rnaseq_sample_metadata_df = get_joined_metadata_df()

    analysis_batches = set([b for b in ANALYSIS_BATCHES.keys() if b])
    star_pipeline_batches = set([b for b in rnaseq_sample_metadata_df["star_pipeline_batch"] if b])

    p = batch_utils.init_arg_parser(default_cpu=1, gsa_key_file=os.path.expanduser("~/.config/gcloud/misc-270914-cb9992ec9b25.json"))
    p.add_argument("--normalize-read-counts", action="store_true", help="whether to normalize unique- and multi-mapped read counts rather than just summing them across input tables")
    p.add_argument("--save-individual-tables", action="store_true", help="Also export individual .bed files with additional columns")
    p.add_argument("batch_name", nargs="+", choices=analysis_batches | star_pipeline_batches, help="Name of RNA-seq batch to process")
    args = p.parse_args()

    if not args.force:
        hl.init(log="/dev/null", quiet=True)

    # process batches
    batch_label = args.batch_name[0] if len(args.batch_name) == 1 else f"{len(args.batch_name)} batches"
    with batch_utils.run_batch(args, batch_name=f"combine junctions: {batch_label}") as batch:
        for batch_name in args.batch_name:
            if batch_name in star_pipeline_batches:
                output_dir = f"gs://macarthurlab-rnaseq/{batch_name}/combined_SJ_out_tables/"
                SJ_out_tab_paths = list(rnaseq_sample_metadata_df[rnaseq_sample_metadata_df["star_pipeline_batch"] == batch_name].star_SJ_out_tab)
            elif batch_name in analysis_batches:
                output_dir = f"gs://macarthurlab-rnaseq/combined_SJ_out_tables/{batch_name}/"
                SJ_out_tab_paths = rnaseq_sample_metadata_df[rnaseq_sample_metadata_df["sample_id"].isin(ANALYSIS_BATCHES[batch_name]["samples"])].star_SJ_out_tab
            else:
                p.error(f"Unexpected batch name: {batch_name}")

            combine_splice_junctions(
                    args,
                    batch,
                    batch_name,
                    SJ_out_tab_paths,
                    args.save_individual_tables,
                    args.normalize_read_counts,
                    output_dir)


    # columns: SAMPID, SMTS (tissue), SMTSD (tissue detail), SMRDLGTH (read length), not stranded, SMRIN (rin), SEX, AGE, DTHHRDY (cause of death)
    #gtex_rnaseq_sample_metadata_df = get_gtex_rnaseq_sample_metadata_df()
