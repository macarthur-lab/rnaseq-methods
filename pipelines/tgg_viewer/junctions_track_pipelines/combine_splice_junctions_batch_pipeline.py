import hail as hl
import logging
import os

from batch import batch_utils
from sample_metadata.rnaseq_metadata_utils import get_rnaseq_downstream_analysis_metadata_df

logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s', level=logging.INFO)
logger = logging.getLogger(__name__)


GCLOUD_PROJECT = "seqr-project"
GCLOUD_USER_ACCOUNT = "weisburd@broadinstitute.org"
GCLOUD_CREDENTIALS_LOCATION = "gs://weisburd-misc/creds"

DOCKER_IMAGE = "weisburd/junctions-track-pipeline@sha256:9f5258ee03c68cea7fd30d8f494eb81aaba883b68e068d6a6e3fe02188bc36a7"

hl.init(log="/dev/null")
batch_utils.set_gcloud_project(GCLOUD_PROJECT)


def combine_splice_junctions(args, batch, batch_name, SJ_out_tab_paths, output_dir):
    normalized_suffix = ".normalized" if args.normalize_read_counts else ""
    gtex_only_prefix = "gtex_only_" if args.gtex_only else ""
    output_filename = f"combined.{batch_name}.{len(SJ_out_tab_paths)}_{gtex_only_prefix}samples{normalized_suffix}.SJ.out.tsv.gz"
    output_path = os.path.join(output_dir, output_filename)
    output_path_exists = hl.hadoop_is_file(output_path)

    output_filename2 = output_filename.replace(".SJ.out.tsv.gz", ".junctions.bed.gz")
    output_path2 = os.path.join(output_dir, output_filename2)
    output_path2_exists = hl.hadoop_is_file(output_path2)

    if not args.force and output_path_exists and output_path2_exists:
        logger.info(f"Output files \n{output_path} and \n{output_path2} exist.\n Skipping...")
        return

    j = batch_utils.init_job(batch, f"combine junctions ({normalized_suffix}): {batch_name} ({len(SJ_out_tab_paths)} files)", DOCKER_IMAGE if not args.raw else None, args.cpu, args.cpu*3.75)
    batch_utils.switch_gcloud_auth_to_user_account(j, GCLOUD_CREDENTIALS_LOCATION, GCLOUD_USER_ACCOUNT, GCLOUD_PROJECT)

    if args.force or not output_path_exists:
        local_SJ_out_tab_paths = []
        for SJ_out_tab_path in SJ_out_tab_paths:
            local_path = batch_utils.localize_file(j, SJ_out_tab_path, use_gcsfuse=False)
            local_SJ_out_tab_paths.append(local_path)

        #local_gencode_gff_path = batch_utils.localize_file(j, "gencode.v26.annotation.gff3.gz", use_gcsfuse=False)

        save_individual_tables_option = "--save-individual-tables" if args.save_individual_tables else ""
        normalize_read_counts_option = "--normalize-read-counts" if args.normalize_read_counts else ""

        local_SJ_out_tab_paths = " ".join(local_SJ_out_tab_paths)
        j.command(f"python3 -u combine_splice_junctions_using_pandas.py --add-sample-id-column -n 20 "
            f"{save_individual_tables_option} "
            f"{normalize_read_counts_option} "
            f"{local_SJ_out_tab_paths}")
        j.command(f"mv combined.{len(SJ_out_tab_paths)}_samples{normalized_suffix}.SJ.out.tsv.gz {output_filename}")
        j.command(f"""gsutil -u {GCLOUD_PROJECT} -m cp *.tsv.gz {output_dir}""")

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
        j.command(f"""gsutil -u {GCLOUD_PROJECT} -m cp {output_filename2}* {output_dir}""")
    else:
        logger.info(f"Output file {output_path2} exists. Skipping...")

    logger.info(f"Output: {output_path}")
    logger.info(f"Output2: {output_path2}")


if __name__ == "__main__":
    # columns: sample_id, star_pipeline_batch, star_SJ_out_tab, 'imputed sex', 'imputed tissue', 'stranded? (rnaseqc)', 'read length (rnaseqc)'
    rnaseq_sample_metadata_df = get_rnaseq_downstream_analysis_metadata_df()
    tissues = set(rnaseq_sample_metadata_df["tissue"])
    batches = set([b for b in rnaseq_sample_metadata_df["batch"] if b])

    p = batch_utils.init_arg_parser(default_cpu=16, gsa_key_file=os.path.expanduser("~/.config/gcloud/misc-270914-cb9992ec9b25.json"))
    p.add_argument("--normalize-read-counts", action="store_true", help="whether to normalize unique- and multi-mapped read counts rather than just summing them across input tables")
    p.add_argument("--save-individual-tables", action="store_true", help="Also export individual .bed files with additional columns")
    g = p.add_mutually_exclusive_group()
    g.add_argument("--include-gtex", action="store_true", help="Include GTEx samples")
    g.add_argument("--gtex-only", action="store_true", help="Include only GTEx samples")
    p.add_argument("batch_name", nargs="+", choices=tissues | batches, help="Name of RNA-seq batch to process")
    args = p.parse_args()

    if not args.force:
        hl.init(log="/dev/null", quiet=True)

    # process batches
    batch_label = args.batch_name[0] if len(args.batch_name) == 1 else f"{len(args.batch_name)} batches ({'normalized' if args.normalize_read_counts else ''})"
    with batch_utils.run_batch(args, batch_name=f"combine junctions: {batch_label}") as batch:
        for batch_name in args.batch_name:
            if batch_name in batches:
                output_dir = f"gs://macarthurlab-rnaseq/{batch_name}/combined_SJ_out_tables/"
                SJ_out_tab_paths = list(rnaseq_sample_metadata_df[rnaseq_sample_metadata_df["batch"] == batch_name].star_SJ_out_tab)
            elif batch_name in tissues:
                output_dir = f"gs://macarthurlab-rnaseq/combined_SJ_out_tables/{batch_name}/"
                SJ_out_tab_paths = rnaseq_sample_metadata_df[rnaseq_sample_metadata_df["tissue"] == batch_name].star_SJ_out_tab
            else:
                p.error(f"Unexpected batch name: {batch_name}")

            if args.gtex_only:
                SJ_out_tab_paths = [p for p in SJ_out_tab_paths if "gtex" in p.lower()]
            elif not args.include_gtex:
                SJ_out_tab_paths = [p for p in SJ_out_tab_paths if "gtex" not in p.lower()]

            combine_splice_junctions(
                    args,
                    batch,
                    batch_name,
                    SJ_out_tab_paths,
                    output_dir)


    # columns: SAMPID, SMTS (tissue), SMTSD (tissue detail), SMRDLGTH (read length), not stranded, SMRIN (rin), SEX, AGE, DTHHRDY (cause of death)
    #gtex_rnaseq_sample_metadata_df = get_gtex_rnaseq_sample_metadata_df()
