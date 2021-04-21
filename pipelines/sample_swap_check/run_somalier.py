import hail as hl  # used for hadoop file utils
import logging
import os
import pandas as pd

from sample_metadata.rnaseq_metadata_utils import get_rnaseq_metadata_joined_with_paths_df

from batch import batch_utils

hl.init(log="/dev/null")

logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s', level=logging.INFO)
logger = logging.getLogger(__name__)

DOCKER_IMAGE = "weisburd/somalier@sha256:4823c2cea6e803443a4d2f31fdff75b315f9fe0934fceee091e1fed2de6ec170"
GCLOUD_PROJECT = "seqr-project"
GCLOUD_USER_ACCOUNT = "weisburd@broadinstitute.org"
GCLOUD_CREDENTIALS_LOCATION = "gs://weisburd-misc/creds"

CPU=0.25

def transfer_metadata_columns_from_df(samples_df, source_df):
    df = pd.DataFrame()
    df.loc[:, 'sample_id'] = source_df['sample_id']
    df.loc[source_df.sample_id, 'batch'] = source_df['star_pipeline_batch']
    df.loc[source_df.sample_id, 'star_bam'] = source_df['star_bam']
    df.loc[source_df.sample_id, 'star_bai'] = source_df['star_bai']
    df.loc[source_df.sample_id, 'grch38_vcf'] = source_df['grch38_vcf']

    df.loc[source_df.sample_id, 'output_dir'] = source_df['star_pipeline_batch'].apply(
        lambda batch_name: f"gs://macarthurlab-rnaseq/{batch_name}/somalier_sample_swap")

    return pd.concat([samples_df, df], axis="rows")


def main():
    rnaseq_sample_metadata_df = get_rnaseq_metadata_joined_with_paths_df()
    #gtex_rnaseq_sample_metadata_df = get_gtex_rnaseq_sample_metadata_df()

    p = batch_utils.init_arg_parser(default_cpu=4, gsa_key_file=os.path.expanduser("~/.config/gcloud/misc-270914-cb9992ec9b25.json"))
    grp = p.add_mutually_exclusive_group(required=True)
    grp.add_argument("--all", action="store_true", help="Process all samples")
    grp.add_argument("-b", "--rnaseq-batch-name", nargs="*", help="RNA-seq batch names to process (eg. -b batch1 batch2)",
        choices=set(rnaseq_sample_metadata_df['star_pipeline_batch']) | set(["gtex_muscle", "gtex_fibroblasts", "gtex_blood"]))
    grp.add_argument("-s", "--rnaseq-sample-id", nargs="*", help="RNA-seq sample IDs to process (eg. -s sample1 sample2)",
        choices=set(rnaseq_sample_metadata_df['sample_id']) | set(['GTEX-1LG7Z-0005-SM-DKPQ6', 'GTEX-PX3G-0006-SM-5SI7E', 'GTEX-1KXAM-0005-SM-DIPEC']))
    args = p.parse_args()

    #rnaseq_sample_metadata_df = rnaseq_sample_metadata_df[
    #    (rnaseq_sample_metadata_df['star_bam'].str.len() > 1) &
    #    (rnaseq_sample_metadata_df['grch38_vcf'].str.len() > 1)
    #]

    # Generate samples_df with these columns: sample_id, bam_path, bai_path, output_dir, batch_name, sex, RIN, ancestry, etc.
    samples_df = pd.DataFrame()
    if args.all:
        samples_df = transfer_metadata_columns_from_df(samples_df, rnaseq_sample_metadata_df)
    elif args.rnaseq_batch_name:
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
            input_bam = metadata_row['star_bam']
            input_vcf = metadata_row['grch38_vcf']
            output_dir = metadata_row['output_dir']

            print("Input bam: ", input_bam)
            output_filename = f"{sample_id}.somalier_results.tsv"
            output_file_path = os.path.join(output_dir, output_filename)

            # check if output file already exists
            if hl.hadoop_is_file(output_file_path) and not args.force:
                logger.info(f"{sample_id} output file already exists: {output_file_path}. Skipping...")
                continue

            j = batch_utils.init_job(batch, f"somalier: {sample_id}", cpu=CPU, memory=CPU*3.75, image=DOCKER_IMAGE)
            batch_utils.switch_gcloud_auth_to_user_account(j, GCLOUD_CREDENTIALS_LOCATION, GCLOUD_USER_ACCOUNT, GCLOUD_PROJECT)

            local_fasta = batch_utils.localize_file(j, batch_utils.HG38_REF_PATHS.fasta, use_gcsfuse=True)
            local_bam_path = batch_utils.localize_file(j, input_bam, use_gcsfuse=True)

            j.command(f"pwd && ls")

            # make sure {sample_id}.groups.tsv exists even if somalier fails to create it
            j.command(f"touch somalier.groups.tsv")
            if input_vcf and len(input_vcf.strip()) > 1:
                local_vcf_path = batch_utils.localize_file(j, input_vcf, use_gcsfuse=True)
                j.command(f"time somalier extract -f {local_fasta} -s sites.hg38.rna.vcf.gz -d output {local_bam_path}")
                j.command(f"time somalier extract -f {local_fasta} -s sites.hg38.rna.vcf.gz -d output {local_vcf_path}")
                j.command(f"time somalier relate output/*.somalier")

                for key in "samples", "pairs", "groups":
                    j.command(f"gsutil -u {GCLOUD_PROJECT} -m cp somalier.{key}.tsv {output_dir}/{sample_id}.{key}.tsv")

            j.command(f"samtools view -c {local_bam_path} chrX > {sample_id}.chrX_read_count.txt")
            j.command(f"samtools view -c {local_bam_path} chrY > {sample_id}.chrY_read_count.txt")
            j.command(f"python3 -c \"print(float(open('{sample_id}.chrX_read_count.txt').read())/float(open('{sample_id}.chrY_read_count.txt').read()))\" > {sample_id}.chrX_over_chrY_read_counts.txt")

            j.command(f"echo -e \"sample_id\\trelatedness\\tchrX/chrY\" >> {sample_id}.somalier_results.tsv")
            j.command(f"echo -e \"{sample_id}\\t$(cut -f 2 somalier.groups.tsv)\\t$(cut -f 1 {sample_id}.chrX_over_chrY_read_counts.txt)\" >> {sample_id}.somalier_results.tsv")
            j.command(f"ls -l")

            j.command(f"gsutil -u {GCLOUD_PROJECT} -m cp *.txt {output_dir}/")
            j.command(f"gsutil -u {GCLOUD_PROJECT} -m cp {sample_id}.somalier_results.tsv {output_dir}/")

            j.command(f"echo Done: {output_file_path}")
            j.command(f"date")

            print("Output file path: ", output_file_path)


if __name__ == "__main__":
    main()
