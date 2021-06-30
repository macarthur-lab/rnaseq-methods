import hail as hl  # used for hadoop file utils
import logging
import os

from batch import batch_utils
from sample_metadata.rnaseq_metadata_utils import get_rnaseq_metadata_joined_with_paths_df

hl.init(log="/dev/null")

logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s', level=logging.INFO)
logger = logging.getLogger(__name__)

GCLOUD_PROJECT = "seqr-project"
GCLOUD_USER_ACCOUNT = "weisburd@broadinstitute.org"
GCLOUD_CREDENTIALS_LOCATION = "gs://weisburd-misc/creds"

DOCKER_IMAGE = "weisburd/junctions-track-pipeline@sha256:93d5a5b17a0031dcfab8a68bc5c2704b19b3f12137ea975d7a8f7e03a701bbc7"


def combine_bigWigs(args, batch, batch_name, bigWig_paths, output_dir):

    output_filename = f"combined.{batch_name}.{len(bigWig_paths)}_samples.bigWig"
    output_path = os.path.join(output_dir, output_filename)
    output_path_exists = hl.hadoop_is_file(output_path)

    if not args.force and output_path_exists:
        logger.info(f"Output file {output_path} exists. Skipping...")
        return

    j = batch_utils.init_job(batch, f"combine bigWigs: {batch_name} ({len(bigWig_paths)} files)", DOCKER_IMAGE if not args.raw else None, args.cpu, args.cpu*3.75)
    batch_utils.switch_gcloud_auth_to_user_account(j, GCLOUD_CREDENTIALS_LOCATION, GCLOUD_USER_ACCOUNT, GCLOUD_PROJECT)

    local_fasta_fai = batch_utils.localize_file(j, batch_utils.HG38_REF_PATHS.fai, use_gcsfuse=False)

    local_bigWig_paths = []
    for bigWig_path in bigWig_paths:
        local_path = batch_utils.localize_file(j, bigWig_path, use_gcsfuse=False)
        local_bigWig_paths.append(local_path)

    local_bigWig_paths = " ".join(local_bigWig_paths)

    j.command(f"""cut -f 1,2 {local_fasta_fai} > {local_fasta_fai}.chromSizes
bigWigMerge {local_bigWig_paths} combined.bedGraph
LC_COLLATE=C sort -k1,1 -k2,2n ./combined.bedGraph > combined.sorted.bedGraph
rm combined.bedGraph
bedGraphToBigWig combined.sorted.bedGraph {local_fasta_fai}.chromSizes {output_filename}

gsutil -m cp {output_filename} {output_dir}""")

    logger.info(f"Output: {output_path}")


if __name__ == "__main__":
    # columns: sample_id, star_pipeline_batch, star_SJ_out_tab, 'imputed sex', 'imputed tissue', 'stranded? (rnaseqc)', 'read length (rnaseqc)'
    rnaseq_sample_metadata_df = get_rnaseq_metadata_joined_with_paths_df()

    analysis_batches = set([b for b in rnaseq_sample_metadata_df["analysis batch"] if b.strip() and b != "x"])
    star_pipeline_batches = set([b for b in rnaseq_sample_metadata_df["star_pipeline_batch"] if b])

    p = batch_utils.init_arg_parser(default_cpu=16, gsa_key_file=os.path.expanduser("~/.config/gcloud/misc-270914-cb9992ec9b25.json"))
    p.add_argument("batch_name", nargs="+", choices=analysis_batches | star_pipeline_batches, help="Name of RNA-seq batch to process")
    args = p.parse_args()

    # process batches
    batch_label = args.batch_name[0] if len(args.batch_name) == 1 else f"{len(args.batch_name)} batches"
    with batch_utils.run_batch(args, batch_name=f"combine junctions: {batch_label}") as batch:
        for batch_name in args.batch_name:
            if batch_name in star_pipeline_batches:
                output_dir = f"gs://tgg-rnaseq/{batch_name}/combined_bigWigs/"
                bigWig_paths = list(rnaseq_sample_metadata_df[rnaseq_sample_metadata_df["star_pipeline_batch"] == batch_name].coverage_bigwig)
            elif batch_name in analysis_batches:
                output_dir = f"gs://tgg-rnaseq/combined_bigWigs/{batch_name}/"
                bigWig_paths = rnaseq_sample_metadata_df[rnaseq_sample_metadata_df["analysis batch"] == batch_name].coverage_bigwig
            else:
                p.error(f"Unexpected batch name: {batch_name}")

            combine_bigWigs(
                args,
                batch,
                batch_name,
                bigWig_paths,
                output_dir)


    # columns: SAMPID, SMTS (tissue), SMTSD (tissue detail), SMRDLGTH (read length), not stranded, SMRIN (rin), SEX, AGE, DTHHRDY (cause of death)
    #gtex_rnaseq_sample_metadata_df = get_gtex_rnaseq_sample_metadata_df()



"""
export NUM_SAMPLES=$(ls -1 GTEX-*.bigWig | wc -l)
echo "Merging: $NUM_SAMPLES samples"

set -x

let boundary1=${NUM_SAMPLES}/2
let boundary2=${boundary1}+1

# split into 2 so bigWigMerge doesn't run out of memory
bigWigMerge $(ls -1 GTEX-*.bigWig | head -n ${boundary1}) temp1.bedGraph
bigWigMerge $(ls -1 GTEX-*.bigWig | tail -n +${boundary2}) temp2.bedGraph

LC_COLLATE=C sort -k1,1 -k2,2n ./temp1.bedGraph > temp1.sorted.bedGraph
LC_COLLATE=C sort -k1,1 -k2,2n ./temp2.bedGraph > temp2.sorted.bedGraph

rm temp1.bedGraph temp2.bedGraph
bedGraphToBigWig temp1.sorted.bedGraph ~/p1/ref/GRCh38/hg38.fa.chromSizes merged1.bigWig
bedGraphToBigWig temp2.sorted.bedGraph ~/p1/ref/GRCh38/hg38.fa.chromSizes merged2.bigWig

bigWigMerge merged1.bigWig merged2.bigWig temp.bedGraph

LC_COLLATE=C sort -k1,1 -k2,2n ./temp.bedGraph > temp.sorted.bedGraph
bedGraphToBigWig temp.sorted.bedGraph ~/p1/ref/GRCh38/hg38.fa.chromSizes merged.${NUM_SAMPLES}_samples.bigWig


if [ -f "merged.${NUM_SAMPLES}_samples.bigWig" ];  then
    rm temp1.bedGraph temp2.bedGraph temp1.sorted.bedGraph temp2.sorted.bedGraph
    rm merged1.bigWig merged2.bigWig temp.bedGraph
    rm temp.sorted.bedGraph
fi

echo Done
"""
