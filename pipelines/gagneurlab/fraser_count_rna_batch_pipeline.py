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

DOCKER_IMAGE = "weisburd/gagneurlab@sha256:6771fc0341c8eff24ea67835bf21f6fbe6eed6c88401842e77ea7d49aac3ca21"
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

    #j.command("((while true; do uptime; sleep 30; done) & )")
    j.command("set -x")

    return j


def get_sample_set_label(sample_ids):
    byte_string = ", ".join(sorted(sample_ids)).encode()
    h = hashlib.md5(byte_string).hexdigest().upper()
    return f"{len(sample_ids)}_samples_{h[:10]}"


def transfer_metadata_columns_from_GTEx_df(samples_df, source_df, batch_name):
    print("Adding")
    print(source_df[['SAMPID', 'SMRIN']])
    print("----------------")

    df = pd.DataFrame()
    df.loc[:, 'sample_id'] = source_df['SAMPID']
    df.loc[source_df.SAMPID, 'batch'] = batch_name
    df.loc[source_df.SAMPID, 'bam_path'] = source_df['rnaseq_bam']
    df.loc[source_df.SAMPID, 'bai_path'] = source_df['rnaseq_bai']
    df.loc[source_df.SAMPID, 'batch_detail'] = source_df['SMNABTCH']
    df.loc[source_df.SAMPID, 'output_dir'] = f"gs://macarthurlab-rnaseq/gtex_v8/fraser_count_rna/"
    df.loc[source_df.SAMPID, 'project'] = "gtex_v8"
    df.loc[source_df.SAMPID, 'sex'] = source_df['SEX']
    df.loc[source_df.SAMPID, 'age'] = source_df['AGE']
    df.loc[source_df.SAMPID, 'cause_of_death'] = source_df['DTHHRDY']
    df.loc[source_df.SAMPID, 'ancestry'] = None
    df.loc[source_df.SAMPID, 'tissue'] = source_df['SMTS']
    df.loc[source_df.SAMPID, 'tissue_detail'] = source_df['SMTSD']
    df.loc[source_df.SAMPID, 'read_length'] = source_df['SMRDLGTH']
    df.loc[source_df.SAMPID, 'stranded'] = "no"
    df.loc[source_df.SAMPID, 'RIN'] = source_df['SMRIN']

    return pd.concat([samples_df, df], axis="rows")


def transfer_metadata_columns_from_df(samples_df, source_df):
    df = pd.DataFrame()
    df.loc[:, 'sample_id'] = source_df['sample_id']
    df.loc[source_df.sample_id, 'batch'] = source_df['star_pipeline_batch']
    df.loc[source_df.sample_id, 'batch_detail'] = source_df['batch_date_from_hg19_bam_header']
    df.loc[source_df.sample_id, 'bam_path'] = source_df['star_bam']
    df.loc[source_df.sample_id, 'bai_path'] = source_df['star_bai']
    df.loc[source_df.sample_id, 'output_dir'] = source_df['star_pipeline_batch'].apply(lambda batch_name: f"gs://macarthurlab-rnaseq/{batch_name}/fraser_count_rna/")
    df.loc[source_df.sample_id, 'project'] = source_df['proj (seqr)']
    df.loc[source_df.sample_id, 'sex'] = source_df['sex']
    df.loc[source_df.sample_id, 'age'] = None
    df.loc[source_df.sample_id, 'cause_of_death'] = None
    df.loc[source_df.sample_id, 'ancestry'] = None
    df.loc[source_df.sample_id, 'tissue'] = None
    df.loc[source_df.sample_id, 'tissue_detail'] = None
    df.loc[source_df.sample_id, 'read_length'] = source_df['read length (rnaseqc)']
    df.loc[source_df.sample_id, 'stranded'] = source_df['stranded? (rnaseqc)']
    df.loc[source_df.sample_id, 'RIN'] = None

    return pd.concat([samples_df, df], axis="rows", sort=True)


def main():
    rnaseq_sample_metadata_df = get_joined_metadata_df()
    gtex_rnaseq_sample_metadata_df = get_gtex_rnaseq_sample_metadata_df()
    #gtex_rnaseq_sample_metadata_df = gtex_rnaseq_sample_metadata_df #.rename({'SAMPID': 'sample_id'}, axis="columns").set_index('sample_id', drop=False)

    p = argparse.ArgumentParser()
    p.add_argument("--batch-billing-project", default="tgg-rare-disease", help="Batch: billing project. Required if submitting to cluster.")
    p.add_argument("--batch-name", help="Batch: (optional) batch name")
    p.add_argument("--batch-temp-bucket", default="macarthurlab-rnaseq", help="Batch: bucket where it stores temp files. "
        "The batch service-account must have Admin permissions for this bucket. These can be added by running "
        "gsutil iam ch serviceAccount:[SERVICE_ACCOUNT_NAME]:objectAdmin gs://[BUCKET_NAME]")
    p.add_argument("-t", "--cpu", type=float, help="Batch: (optional) number of CPUs (eg. 0.5)", default=1, choices=[0.25, 0.5, 1, 2, 4, 8, 16])
    p.add_argument("-m1", "--memory-step1", type=float, help="Batch: (optional) memory in gigabytes (eg. 3.75)", default=3.75)
    p.add_argument("-m2", "--memory-step2", type=float, help="Batch: (optional) memory in gigabytes (eg. 3.75)", default=3.75)
    p.add_argument("-f", "--force", action="store_true", help="Recompute and overwrite cached or previously computed data")
    p.add_argument("--skip-step1", action="store_true", help="Skip count-split-reads step")
    grp = p.add_mutually_exclusive_group(required=True)
    grp.add_argument("-b", "--rnaseq-batch-name", nargs="*", help="RNA-seq batch names to process (eg. -b batch1 batch2)",
        choices=set(rnaseq_sample_metadata_df['star_pipeline_batch']) | set(["muscle", "fibroblasts", "whole_blood", "lymphocytes"])
    )
    grp.add_argument("-s", "--rnaseq-sample-id", nargs="*", help="RNA-seq sample IDs to process (eg. -s sample1 sample2)",
        choices=set(rnaseq_sample_metadata_df['sample_id']) | set(['GTEX-1LG7Z-0005-SM-DKPQ6', 'GTEX-PX3G-0006-SM-5SI7E', 'GTEX-1KXAM-0005-SM-DIPEC']))
    p.add_argument("-tsv", "--only-generate-tsv", action="store_true", help="Exit after generating metadata tsv")
    args = p.parse_args()

    #logger.info("\n".join(df.columns))
    # Generate samples_df with these columns: sample_id, bam_path, bai_path, output_dir, batch_name, sex, RIN, ancestry, etc.

    samples_df = pd.DataFrame()
    if args.rnaseq_batch_name:
        for batch_name in args.rnaseq_batch_name:
            if batch_name in set(["muscle", "fibroblasts", "whole_blood", "lymphocytes"]):
                if batch_name == "muscle":
                    gtex_df = gtex_rnaseq_sample_metadata_df[gtex_rnaseq_sample_metadata_df.SMTSD == "Muscle - Skeletal"]
                elif batch_name == "fibroblasts":
                    gtex_df = gtex_rnaseq_sample_metadata_df[gtex_rnaseq_sample_metadata_df.SMTSD == "Cells - Cultured fibroblasts"]
                elif batch_name == "lymphocytes":
                    gtex_df = gtex_rnaseq_sample_metadata_df[gtex_rnaseq_sample_metadata_df.SMTSD == "Cells - EBV-transformed lymphocytes"]
                elif batch_name == "whole_blood":
                    gtex_df = gtex_rnaseq_sample_metadata_df[gtex_rnaseq_sample_metadata_df.SMTSD == "Whole Blood"]
                else:
                    p.error(f"Unexpected batch name: {batch_name}")

                gtex_df = gtex_df.sort_values(by='SMRIN', ascending=False)
                samples_df = transfer_metadata_columns_from_GTEx_df(samples_df, gtex_df[:100], batch_name)

                tgg_df = rnaseq_sample_metadata_df[rnaseq_sample_metadata_df['imputed tissue'] == batch_name]
                samples_df = transfer_metadata_columns_from_df(samples_df, tgg_df)
            else:
                df = rnaseq_sample_metadata_df[rnaseq_sample_metadata_df['star_pipeline_batch'] == batch_name]
                samples_df = transfer_metadata_columns_from_df(samples_df, df)

    elif args.rnaseq_sample_id:
        df = rnaseq_sample_metadata_df[rnaseq_sample_metadata_df.sample_id.isin(set(args.rnaseq_sample_id))]
        samples_df = transfer_metadata_columns_from_df(samples_df, df)
    else:
        p.error("Must specify -b or -s")

    sample_set_label = get_sample_set_label(samples_df.sample_id)

    tsv_output_path = f"metadata_{sample_set_label}.tsv"
    samples_df.to_csv(tsv_output_path, sep="\t", index=False)
    print(f"Wrote {len(samples_df)} samples to {tsv_output_path}")

    if args.only_generate_tsv:
        sys.exit(0)



    logger.info(f"Processing {len(samples_df)} sample ids: {', '.join(samples_df.sample_id[:20])}")

    # see https://hail.zulipchat.com/#narrow/stream/223457-Batch-support/topic/auth.20as.20user.20account for more details
    backend = hb.ServiceBackend(billing_project=args.batch_billing_project, bucket=args.batch_temp_bucket)
    b = hb.Batch(backend=backend, name=args.batch_name)

    split_reads_samples = []

    split_reads_output_files = set()
    split_reads_jobs = {}

    non_split_reads_output_files = set()
    non_split_reads_jobs = {}

    j_extract_splice_junctions = None
    #j_calculate_psi_values = None
    for step in 1, 2:
        for sample_id in samples_df.sample_id:
            metadata_row = samples_df.loc[sample_id]

            # set job inputs & outputs
            input_bam, input_bai = metadata_row['bam_path'], metadata_row['bai_path']
            output_dir = metadata_row['output_dir']

            #input_bam = "gs://macarthurlab-rnaseq/temp/MUN_FAM5_SIBLINGMDC1A_01_R1.Aligned.sortedByCoord.out.subset.bam"
            #input_bai = "gs://macarthurlab-rnaseq/temp/MUN_FAM5_SIBLINGMDC1A_01_R1.Aligned.sortedByCoord.out.subset.bam.bai"

            print("Input bam: ", input_bam)
            if step == 1:
                output_file_path = os.path.join(output_dir, f"fraser_count_split_reads_{sample_id}.tar.gz")
                memory = args.memory_step1
            elif step == 2:
                output_file_path = os.path.join(output_dir, f"fraser_count_non_split_reads_{sample_id}__{sample_set_label}.tar.gz")
                memory = args.memory_step2

            if step == 1:
                split_reads_samples.append(sample_id)
                split_reads_output_files.add(output_file_path.replace(sample_id, "*"))
            elif step == 2:
                non_split_reads_output_files.add(output_file_path.replace(sample_id, "*"))

            if step == 1 and args.skip_step1:
                continue

            # check if output file already exists
            if hl.hadoop_is_file(output_file_path) and not args.force:
                logger.info(f"{sample_id} output file already exists: {output_file_path}. Skipping...")
                continue

            file_stats = hl.hadoop_stat(metadata_row['bam_path'])
            bam_size = int(round(file_stats['size_bytes']/10.**9))
            disk_size = bam_size * 2

            job_label = f"Count {'split' if step == 1 else 'non-split'} reads"
            j = init_job(b, name=f"{job_label}: {sample_id}", cpu=args.cpu, memory=memory, disk_size=disk_size, switch_to_user_account=True, image=DOCKER_IMAGE)

            j.command(f"gsutil -u {GCLOUD_PROJECT} -m cp {input_bam} {sample_id}.bam")
            j.command(f"gsutil -u {GCLOUD_PROJECT} -m cp {input_bai} {sample_id}.bam.bai")
            j.command(f"touch {sample_id}.bam.bai")
            bam_path = f"{sample_id}.bam"

            j.command(f"pwd && ls && date")

            if step == 1:
                j.command(f"Rscript --vanilla /countSplitReads.R {sample_id} {bam_path}")
            elif step == 2:
                if sample_id in split_reads_jobs:
                    j.depends_on(split_reads_jobs[sample_id])
                if j_extract_splice_junctions:
                    j.depends_on(j_extract_splice_junctions)

                j.command(f"gsutil -m cp {output_file_path_splice_junctions_RDS} .")
                j.command(f"Rscript --vanilla /countNonSplitReads.R {sample_id} {bam_path} {os.path.basename(output_file_path_splice_junctions_RDS)}")

            j.command(f"ls .")
            j.command(f"tar czf {j.output_tar_gz} cache")

            #j.command(f"ls -lh .")
            j.command(f"echo Done: {output_file_path}")
            j.command(f"date")

            # copy output
            b.write_output(j.output_tar_gz, output_file_path)
            print("Output file path: ", output_file_path)

            if step == 1:
                split_reads_jobs[sample_id] = j
            elif step == 2:
                non_split_reads_jobs[sample_id] = j

        if len(split_reads_output_files) == 0:
            break

        if step == 1:
            output_file_path_splice_junctions_RDS = os.path.join("gs://macarthurlab-rnaseq/fraser/", f"spliceJunctions_{sample_set_label}.RDS")
            if hl.hadoop_is_file(output_file_path_splice_junctions_RDS) and not args.force:
                logger.info(f"{output_file_path_splice_junctions_RDS} file already exists. Skipping extractSpliceJunctions.R step...")
                continue

            j_extract_splice_junctions = init_job(b, name=f"Extract splice-junctions", disk_size=30, memory=64, switch_to_user_account=True, image=DOCKER_IMAGE)
            for j in split_reads_jobs.values():
                j_extract_splice_junctions.depends_on(j)

            j_extract_splice_junctions.command(f"gsutil -m cp {' '.join(split_reads_output_files)} .")
            j_extract_splice_junctions.command(f"gsutil -m cp gs://macarthurlab-rnaseq/fraser/bam_header.bam .")
            j_extract_splice_junctions.command(f"for i in fraser_count_split_reads*.tar.gz; do tar xzf $i; done")
            j_extract_splice_junctions.command(f"pwd && ls && date")
            j_extract_splice_junctions.command(f"Rscript --vanilla /extractSpliceJunctions.R --num-threads=4 bam_header.bam")
            j_extract_splice_junctions.command(f"ls .")
            j_extract_splice_junctions.command(f"gsutil -m cp spliceJunctions.RDS {output_file_path_splice_junctions_RDS}")
            print("Output file path: ", output_file_path_splice_junctions_RDS)
        elif step == 2:
            output_file_path = os.path.join("gs://macarthurlab-rnaseq/fraser/", f"calculatedPSIValues_{sample_set_label}.tar.gz")
            if hl.hadoop_is_file(output_file_path) and not args.force:
                logger.info(f"{output_file_path} file already exists. Skipping calculatePSIValues.R step...")
                continue

            j_calculate_psi_values = init_job(b, name=f"Calculate PSI values", disk_size=50, cpu=16, memory=64, switch_to_user_account=True, image=DOCKER_IMAGE)
            if j_extract_splice_junctions:
                j_calculate_psi_values.depends_on(j_extract_splice_junctions)
            for j in non_split_reads_jobs.values():
                j_calculate_psi_values.depends_on(j)

            j_calculate_psi_values.command(f"gsutil -m cp {' '.join(split_reads_output_files)} .")
            j_calculate_psi_values.command(f"gsutil -m cp {' '.join(non_split_reads_output_files)} .")
            j_calculate_psi_values.command(f"gsutil -m cp {output_file_path_splice_junctions_RDS} .")
            j_calculate_psi_values.command(f"gsutil -m cp gs://macarthurlab-rnaseq/fraser/bam_header.bam .")
            j_calculate_psi_values.command(f"for i in fraser_count_split_reads*.tar.gz; do tar xzf $i; done")
            j_calculate_psi_values.command(f"for i in fraser_count_non_split_reads*.tar.gz; do tar xzf $i; done")
            j_calculate_psi_values.command(f"pwd && ls && date")
            j_calculate_psi_values.command(f"Rscript --vanilla /calculatePSIValues.R --num-threads=4 {os.path.basename(output_file_path_splice_junctions_RDS)} bam_header.bam")
            j_calculate_psi_values.command(f"ls .")
            #j_calculate_psi_values.command(f"cp fdsWithPSIValues.RDS {j_calculate_psi_values.fdsWithPSIValues}")
            j_calculate_psi_values.command(f"tar czf {j_calculate_psi_values.output_tar_gz} cache savedObjects fdsWithPSIValues.RDS")
            b.write_output(j_calculate_psi_values.output_tar_gz, output_file_path)
            print("Output file path: ", output_file_path)

    b.run()

    if isinstance(backend, hb.ServiceBackend):
        backend.close()


if __name__ == "__main__":
    main()
