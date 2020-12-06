import datetime
import hail as hl
import hashlib
import logging
import os
import pandas as pd
import sys


from batch import batch_utils
from sample_metadata.rnaseq_metadata_utils import ANALYSIS_BATCHES
from gagneurlab.gagneur_utils import ALL_METADATA_TSV, BAM_HEADER_PATH, GENCODE_TXDB, DOCKER_IMAGE, GCLOUD_PROJECT, GCLOUD_CREDENTIALS_LOCATION, GCLOUD_USER_ACCOUNT
from gagneurlab.fraser_batch_pipeline_Rscripts import get_EXTRACT_SPLICE_JUNCTIONS_Rscript, \
    get_CALCULATE_PSI_VALUES_Rscript, get_CALCULATE_BEST_Q_Rscript, get_RUN_FRASER_ANALYSIS_Rscript, \
    get_FILTER_AND_ANNOTATE_DATA_Rscript, get_PLOT_RESULTS

logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s', level=logging.INFO)
logger = logging.getLogger(__name__)

# NOTE: xvfb-run is used in the Rscript commands below to solve the issue of R graphics libraries requiring a
# computer monitor and related graphics libraries to be set up in the execution environment

PADJ_THRESHOLD=0.1
DELTA_PSI_THRESHOLD=0.1


def extract_splice_junctions(j_extract_splice_junctions, split_reads_files, num_cpu, output_file_path_splice_junctions_RDS):
    batch_utils.switch_gcloud_auth_to_user_account(j_extract_splice_junctions, GCLOUD_CREDENTIALS_LOCATION, GCLOUD_USER_ACCOUNT, GCLOUD_PROJECT)

    for split_reads_output_files_batch in [split_reads_files[i:i+10] for i in range(0, len(split_reads_files), 10)]:
        j_extract_splice_junctions.command(f"gsutil -m cp {' '.join(split_reads_output_files_batch)} .")
    j_extract_splice_junctions.command(f"gsutil -m cp {BAM_HEADER_PATH} .")
    j_extract_splice_junctions.command(f"for i in fraser_count_split_reads*.tar.gz; do tar xzf $i; done")
    j_extract_splice_junctions.command(f"pwd && ls -lh && date && echo ------- && find cache -name '*.*'")
    j_extract_splice_junctions.command(f"""time xvfb-run Rscript -e '{get_EXTRACT_SPLICE_JUNCTIONS_Rscript(BAM_HEADER_PATH, num_cpu)}'""")
    j_extract_splice_junctions.command(f"ls -lh .")
    j_extract_splice_junctions.command(f"echo ===============; echo cache dir; echo ===============; find cache")
    j_extract_splice_junctions.command(f"echo ===============; echo savedObects dir; echo ===============; find savedObjects")

    j_extract_splice_junctions.command(f"gsutil -m cp spliceJunctions.RDS {output_file_path_splice_junctions_RDS}")

    print("Output file path: ", output_file_path_splice_junctions_RDS)


def calculate_psi_values(j_calculate_psi_values, sample_set_label, split_reads_files, non_split_reads_files, splice_junctions_RDS_path, metadata_tsv_path, output_file_path_calculated_psi_values_tar_gz):
    batch_utils.switch_gcloud_auth_to_user_account(j_calculate_psi_values, GCLOUD_CREDENTIALS_LOCATION, GCLOUD_USER_ACCOUNT, GCLOUD_PROJECT)

    j_calculate_psi_values.command(f"mkdir -p /tmp/fraser/{sample_set_label}")  # work-around for https://github.com/c-mertes/FRASER/issues/11
    j_calculate_psi_values.command(f"cd /tmp/fraser/{sample_set_label}")
    for split_reads_output_files_batch in [split_reads_files[i:i+10] for i in range(0, len(split_reads_files), 10)]:
        j_calculate_psi_values.command(f"gsutil -m cp {' '.join(split_reads_output_files_batch)} .")
    for non_split_reads_output_files_batch in [non_split_reads_files[i:i+10] for i in range(0, len(non_split_reads_files), 10)]:
        j_calculate_psi_values.command(f"gsutil -m cp {' '.join(non_split_reads_output_files_batch)} .")
    j_calculate_psi_values.command(f"gsutil -m cp {splice_junctions_RDS_path} .")
    j_calculate_psi_values.command(f"gsutil -m cp {metadata_tsv_path} .")

    j_calculate_psi_values.command(f"gsutil -m cp {BAM_HEADER_PATH} .")

    j_calculate_psi_values.command(f"for i in fraser_count_split_reads*.tar.gz; do tar xzf $i; done")
    j_calculate_psi_values.command(f"for i in fraser_count_non_split_reads*.tar.gz; do tar xzf $i; done")
    j_calculate_psi_values.command(f"rm cache/nonSplicedCounts/Data_Analysis/spliceSiteCoordinates.RDS")
    j_calculate_psi_values.command(f"pwd && ls -lh && date && echo ------- && find cache -name '*.*'")
    j_calculate_psi_values.command(f"""time xvfb-run Rscript -e '{get_CALCULATE_PSI_VALUES_Rscript(splice_junctions_RDS_path, metadata_tsv_path, BAM_HEADER_PATH, num_cpu=1)}'""")
    j_calculate_psi_values.command(f"pwd")
    j_calculate_psi_values.command(f"ls -lh .")
    j_calculate_psi_values.command(f"echo ===============; echo cache dir; echo ===============; find cache")
    j_calculate_psi_values.command(f"echo ===============; echo savedObects dir; echo ===============; find savedObjects")
    j_calculate_psi_values.command(f"rm *.tar.gz *.bam")
    j_calculate_psi_values.command(f"cd ..")
    j_calculate_psi_values.command(f"tar czf {os.path.basename(output_file_path_calculated_psi_values_tar_gz)} {sample_set_label}")
    j_calculate_psi_values.command(f"gsutil -m cp {os.path.basename(output_file_path_calculated_psi_values_tar_gz)} {output_file_path_calculated_psi_values_tar_gz}")
    print("Output file path: ", output_file_path_calculated_psi_values_tar_gz)


def filter_and_annotate_data(j_filter_and_annotate_data, sample_set_label, calculated_psi_values_tar_gz_path, output_file_path_filter_and_annotate_data_tar_gz):
    batch_utils.switch_gcloud_auth_to_user_account(j_filter_and_annotate_data, GCLOUD_CREDENTIALS_LOCATION, GCLOUD_USER_ACCOUNT, GCLOUD_PROJECT)

    j_filter_and_annotate_data.command(f"mkdir -p /tmp/fraser/")  # work-around for https://github.com/c-mertes/FRASER/issues/11
    j_filter_and_annotate_data.command(f"cd /tmp/fraser/")
    j_filter_and_annotate_data.command(f"gsutil -m cp {calculated_psi_values_tar_gz_path} .")
    j_filter_and_annotate_data.command(f"tar xzf {os.path.basename(calculated_psi_values_tar_gz_path)}")
    j_filter_and_annotate_data.command(f"cd {sample_set_label}")
    j_filter_and_annotate_data.command(f"pwd && ls -lh && date && echo ------- && find cache -name '*.*'")
    j_filter_and_annotate_data.command(f"""time xvfb-run Rscript -e '{get_FILTER_AND_ANNOTATE_DATA_Rscript(sample_set_label, delta_psi_threshold=DELTA_PSI_THRESHOLD)}'""")
    j_filter_and_annotate_data.command(f"echo ===============; echo ls .; echo ===============; pwd; ls -lh .")
    j_filter_and_annotate_data.command(f"echo ===============; echo cache dir; echo ===============; find cache")
    j_filter_and_annotate_data.command(f"echo ===============; echo savedObects dir; echo ===============; find savedObjects")
    j_filter_and_annotate_data.command(f"cd ..")
    j_filter_and_annotate_data.command(f"tar czf {os.path.basename(output_file_path_filter_and_annotate_data_tar_gz)} {sample_set_label}")
    j_filter_and_annotate_data.command(f"gsutil -m cp {os.path.basename(output_file_path_filter_and_annotate_data_tar_gz)} {output_file_path_filter_and_annotate_data_tar_gz}")
    print("Output file path: ", output_file_path_filter_and_annotate_data_tar_gz)


def calculate_best_q(j_calculate_best_q, sample_set_label, num_cpu, filter_and_annotate_data_tar_gz_path, output_file_path_calculated_best_q_tar_gz):
    batch_utils.switch_gcloud_auth_to_user_account(j_calculate_best_q, GCLOUD_CREDENTIALS_LOCATION, GCLOUD_USER_ACCOUNT, GCLOUD_PROJECT)

    j_calculate_best_q.command(f"mkdir -p /tmp/fraser/")  # work-around for https://github.com/c-mertes/FRASER/issues/11
    j_calculate_best_q.command(f"cd /tmp/fraser/")
    j_calculate_best_q.command(f"gsutil -m cp {filter_and_annotate_data_tar_gz_path} .")
    j_calculate_best_q.command(f"tar xzf {os.path.basename(filter_and_annotate_data_tar_gz_path)}")
    j_calculate_best_q.command(f"cd {sample_set_label}")
    j_calculate_best_q.command(f"pwd && ls -lh && date && echo ------- && find cache -name '*.*'")
    j_calculate_best_q.command(f"""time xvfb-run Rscript -e '{get_CALCULATE_BEST_Q_Rscript(sample_set_label, num_cpu)}'""")
    j_calculate_best_q.command(f"echo ===============; echo ls .; echo ===============; pwd; ls -lh .")
    j_calculate_best_q.command(f"echo ===============; echo cache dir; echo ===============; find cache")
    j_calculate_best_q.command(f"echo ===============; echo savedObects dir; echo ===============; find savedObjects")
    j_calculate_best_q.command(f"cd ..")
    j_calculate_best_q.command(f"tar czf {os.path.basename(output_file_path_calculated_best_q_tar_gz)} {sample_set_label}")
    j_calculate_best_q.command(f"gsutil -m cp {os.path.basename(output_file_path_calculated_best_q_tar_gz)} {output_file_path_calculated_best_q_tar_gz}")
    print("Output file path: ", output_file_path_calculated_best_q_tar_gz)


def run_fraser_analysis(j_run_fraser_analysis, sample_set_label, calculated_best_q_tar_gz_path, output_file_path_fraser_analysis_tar_gz, output_file_path_fraser_analysis_results_only_tar_gz):
    batch_utils.switch_gcloud_auth_to_user_account(j_run_fraser_analysis, GCLOUD_CREDENTIALS_LOCATION, GCLOUD_USER_ACCOUNT, GCLOUD_PROJECT)

    j_run_fraser_analysis.command(f"mkdir -p /tmp/fraser/")  # work-around for https://github.com/c-mertes/FRASER/issues/11
    j_run_fraser_analysis.command(f"cd /tmp/fraser/")
    j_run_fraser_analysis.command(f"gsutil -m cp {calculated_best_q_tar_gz_path} .")
    j_run_fraser_analysis.command(f"tar xzf {os.path.basename(calculated_best_q_tar_gz_path)}")
    j_run_fraser_analysis.command(f"cd {sample_set_label}")
    j_run_fraser_analysis.command(f"pwd && ls -lh && date && echo ------- && find cache -name '*.*'")
    j_run_fraser_analysis.command(f"""time xvfb-run Rscript -e '{get_RUN_FRASER_ANALYSIS_Rscript(sample_set_label, num_cpu=1, delta_psi_threshold=DELTA_PSI_THRESHOLD, padj_threshold=PADJ_THRESHOLD)}'""")
    j_run_fraser_analysis.command(f"echo ===============; echo ls .; echo ===============; pwd; ls -lh .")
    j_run_fraser_analysis.command(f"cd ..")
    j_run_fraser_analysis.command(f"tar czf {os.path.basename(output_file_path_fraser_analysis_tar_gz)} {sample_set_label}")
    j_run_fraser_analysis.command(f"gsutil -m cp {os.path.basename(output_file_path_fraser_analysis_tar_gz)} {output_file_path_fraser_analysis_tar_gz}")
    print("Output file path: ", output_file_path_fraser_analysis_tar_gz)

    j_run_fraser_analysis.command(f"rm -rf {sample_set_label}/cache {sample_set_label}/savedObjects")
    j_run_fraser_analysis.command(f"tar czf {os.path.basename(output_file_path_fraser_analysis_results_only_tar_gz)} {sample_set_label}")
    j_run_fraser_analysis.command(f"gsutil -m cp {os.path.basename(output_file_path_fraser_analysis_results_only_tar_gz)} {output_file_path_fraser_analysis_results_only_tar_gz}")
    print("Output file path (results only): ", output_file_path_fraser_analysis_results_only_tar_gz)


def run_plots(j_run_fraser_analysis, sample_set_label, fraser_analysis_results_only_tar_gz_path, output_file_path_fraser_volcano_plots_tar_gz):
    batch_utils.switch_gcloud_auth_to_user_account(j_run_fraser_analysis, GCLOUD_CREDENTIALS_LOCATION, GCLOUD_USER_ACCOUNT, GCLOUD_PROJECT)

    j_run_fraser_analysis.command(f"mkdir -p /tmp/fraser/")  # work-around for https://github.com/c-mertes/FRASER/issues/11
    j_run_fraser_analysis.command(f"cd /tmp/fraser/")
    j_run_fraser_analysis.command(f"gsutil -m cp {fraser_analysis_results_only_tar_gz_path} .")
    j_run_fraser_analysis.command(f"tar xzf {os.path.basename(fraser_analysis_results_only_tar_gz_path)}")
    j_run_fraser_analysis.command(f"cd {sample_set_label}")
    j_run_fraser_analysis.command(f"pwd && ls -lh && date && echo ------- && find cache -name '*.*'")
    j_run_fraser_analysis.command(f"""time xvfb-run Rscript -e '{get_PLOT_RESULTS(sample_set_label, delta_psi_threshold=DELTA_PSI_THRESHOLD, padj_threshold=PADJ_THRESHOLD)}'""")
    j_run_fraser_analysis.command(f"echo ===============; echo ls .; echo ===============; pwd; ls -lh .")
    j_run_fraser_analysis.command(f"cd ..")
    j_run_fraser_analysis.command(f"tar czf {os.path.basename(output_file_path_fraser_volcano_plots_tar_gz)} {sample_set_label}")
    j_run_fraser_analysis.command(f"gsutil -m cp {os.path.basename(output_file_path_fraser_volcano_plots_tar_gz)} {output_file_path_fraser_volcano_plots_tar_gz}")
    print("Output file path: ", output_file_path_fraser_volcano_plots_tar_gz)


def main():
    p = batch_utils.init_arg_parser(default_cpu=4, gsa_key_file=os.path.expanduser("~/.config/gcloud/misc-270914-cb9992ec9b25.json"))
    p.add_argument("--with-gtex", help="Use GTEX controls.", action="store_true")
    p.add_argument("--skip-step1", action="store_true", help="Skip count-split-reads step")
    p.add_argument("--skip-step2", action="store_true", help="Skip compute-PSI step")
    p.add_argument("--skip-step3", action="store_true", help="Skip filter-and-annotate step")
    p.add_argument("--skip-step4", action="store_true", help="Skip compute-best-Q step")
    p.add_argument("--skip-step5", action="store_true", help="Skip FRASER analysis step")
    p.add_argument("-m1", "--memory-step1", type=float, help="Batch: (optional) memory in gigabytes (eg. 3.75)", default=3.75)
    p.add_argument("-m2", "--memory-step2", type=float, help="Batch: (optional) memory in gigabytes (eg. 3.75)", default=3.75)
    p.add_argument("--metadata-tsv-path", default=ALL_METADATA_TSV, help="Table with columns: sample_id, bam_path, bai_path, batch")
    p.add_argument("batch_name", nargs="+", choices=ANALYSIS_BATCHES.keys(), help="Name of RNA-seq batch to process")
    args = p.parse_args()

    hl.init(log="/dev/null", quiet=True)

    local_metadata_tsv_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.path.basename(ALL_METADATA_TSV))
    print(f"Copying {local_metadata_tsv_path} to {ALL_METADATA_TSV}")
    hl.hadoop_copy(local_metadata_tsv_path, ALL_METADATA_TSV)
    metadata_tsv_df = pd.read_table(local_metadata_tsv_path)


    with hl.hadoop_open(args.metadata_tsv_path) as f:
        samples_df_unmodified = pd.read_table(f).set_index("sample_id", drop=False)

    batch_label = f"FRASER"
    if args.with_gtex:
        batch_label += f" (with GTEx)"
    batch_label += ": "
    batch_label += ','.join(args.batch_name)
    batch_label += f" (dPsi={DELTA_PSI_THRESHOLD}, padj={PADJ_THRESHOLD})"
    with batch_utils.run_batch(args, batch_label) as batch:

        for batch_name in args.batch_name:
            samples_df = samples_df_unmodified
            batch_dict = ANALYSIS_BATCHES[batch_name]
            batch_tissue = batch_dict['tissue']
            batch_sex = batch_dict['sex']

            num_sample_missing_from_metadata_tsv = len(batch_dict['samples']) - len(metadata_tsv_df.sample_id.isin(batch_dict['samples']))
            if num_sample_missing_from_metadata_tsv > 0:
                raise ValueError(f"{batch_tissue} metadata_tsv_df is missing {num_sample_missing_from_metadata_tsv} sample ids. Download count files and run export_gagneur_metadata_table.py")

            sample_ids = list(batch_dict['samples'])
            if args.with_gtex:
                batch_name += "_with_GTEX"
                samples_df_filter = (samples_df.tissue == batch_tissue)
                samples_df_filter &= samples_df.sample_id.str.startswith("GTEX")
                if batch_sex == "M" or batch_sex == "F":
                    samples_df_filter &= (samples_df.sex == batch_sex)
                sample_ids += list(samples_df[samples_df_filter].sample_id)
            else:
                batch_name += "_without_GTEX"

            samples_df = samples_df.loc[sample_ids]
            byte_string = ", ".join(sorted(samples_df.sample_id)).encode()
            h = hashlib.md5(byte_string).hexdigest().upper()
            sample_set_label = f"{batch_name}__{len(samples_df.sample_id)}_samples_{h[:10]}"

            logger.info(f"Processing {sample_set_label}: {len(samples_df)} sample ids: {', '.join(samples_df.sample_id[:20])}")

            split_reads_samples = []

            split_reads_output_files = []
            split_reads_jobs = {}

            non_split_reads_output_files = []
            non_split_reads_jobs = {}

            j_extract_splice_junctions = None
            j_calculate_psi_values = None
            j_filter_and_annotate_data = None
            j_calculate_best_q = None
            j_fraser_analysis = None
            j_get_plots = None

            # based on docs @ https://bioconductor.org/packages/devel/bioc/vignettes/FRASER/inst/doc/FRASER.pdf
            # step 1: count spliced reads
            # step 2: count non-spliced reads at acceptors & donors of splice junctions detected in step 1
            for step in 1, 2:
                for sample_id in samples_df.sample_id:
                    metadata_row = samples_df.loc[sample_id]

                    # set job inputs & outputs
                    input_bam, input_bai = metadata_row['bam_path'], metadata_row['bai_path']
                    if "GTEX" in sample_id:
                        output_dir_for_sample_specific_data = "gs://macarthurlab-rnaseq/gtex_v8/fraser_count_rna/"
                    else:
                        output_dir_for_sample_specific_data = f"gs://macarthurlab-rnaseq/{metadata_row['batch']}/fraser_count_rna/"

                    output_dir_for_batch_specific_data = f"gs://macarthurlab-rnaseq/gagneur/fraser/results/{sample_set_label}"

                    output_file_path_splice_junctions_RDS = os.path.join(output_dir_for_batch_specific_data, f"spliceJunctions_{sample_set_label}.RDS")
                    output_file_path_calculated_psi_values_tar_gz = os.path.join(output_dir_for_batch_specific_data, f"calculatedPSIValues_{sample_set_label}.tar.gz")
                    output_file_path_filter_and_annotate_data_tar_gz = os.path.join(output_dir_for_batch_specific_data, f"filteredAndAnnotated_{sample_set_label}__dpsi_{DELTA_PSI_THRESHOLD}.tar.gz")
                    output_file_path_calculated_best_q_tar_gz = os.path.join(output_dir_for_batch_specific_data, f"calculatedBestQ_{sample_set_label}__dpsi_{DELTA_PSI_THRESHOLD}.tar.gz")
                    output_file_path_fraser_analysis_tar_gz = os.path.join(output_dir_for_batch_specific_data, f"fraserAnalysis_using_PCA_{sample_set_label}__dpsi_{DELTA_PSI_THRESHOLD}_padj_{PADJ_THRESHOLD}.tar.gz")
                    output_file_path_fraser_analysis_results_only_tar_gz = os.path.join(output_dir_for_batch_specific_data, f"fraserAnalysis_using_PCA_{sample_set_label}_results_only__dpsi_{DELTA_PSI_THRESHOLD}_padj_{PADJ_THRESHOLD}.tar.gz")
                    output_file_path_fraser_volcano_plots_tar_gz = os.path.join(output_dir_for_batch_specific_data, f"fraserVolcanoPlots_{sample_set_label}__dpsi_{DELTA_PSI_THRESHOLD}_padj_{PADJ_THRESHOLD}.tar.gz")

                    print("Input bam: ", input_bam)
                    if step == 1:
                        output_file_path = os.path.join(output_dir_for_sample_specific_data, f"fraser_count_split_reads_{sample_id}.tar.gz")
                        memory = args.memory_step1
                    elif step == 2:
                        output_file_path = os.path.join(output_dir_for_batch_specific_data, f"fraser_count_non_split_reads_{sample_id}__{sample_set_label}.tar.gz")
                        memory = args.memory_step2

                    if step == 1:
                        split_reads_samples.append(sample_id)
                        split_reads_output_files.append(output_file_path)
                    elif step == 2:
                        non_split_reads_output_files.append(output_file_path)

                    if (step == 1 and args.skip_step1) or (step == 2 and args.skip_step2):
                        continue

                    # check if output file already exists
                    if not args.force and hl.hadoop_is_file(output_file_path):
                        logger.info(f"{sample_id} output file already exists: {output_file_path}. Skipping...")
                        continue

                    if not args.local:
                        file_stats = hl.hadoop_stat(metadata_row['bam_path'])
                        bam_size = int(round(file_stats['size_bytes']/10.**9))
                        disk_size = bam_size * 2
                    else:
                        disk_size = None

                    job_label = f"Count {'split' if step == 1 else 'non-split'} reads"
                    j = batch_utils.init_job(batch, f"{job_label}: {sample_id}", cpu=args.cpu, memory=memory, disk_size=disk_size, image=DOCKER_IMAGE)
                    batch_utils.switch_gcloud_auth_to_user_account(j, GCLOUD_CREDENTIALS_LOCATION, GCLOUD_USER_ACCOUNT, GCLOUD_PROJECT)

                    j.command(f"gsutil -u {GCLOUD_PROJECT} -m cp {input_bam} {sample_id}.bam")
                    j.command(f"gsutil -u {GCLOUD_PROJECT} -m cp {input_bai} {sample_id}.bam.bai")
                    j.command(f"touch {sample_id}.bam.bai")
                    bam_path = f"{sample_id}.bam"

                    j.command(f"pwd && ls -lh && date")

                    if step == 1:
                        # count split reads
                        j.command(f"""time xvfb-run Rscript -e '
library(FRASER)
library(data.table)

sampleTable = data.table(sampleID=c("{sample_id}"), bamFile=c("{bam_path}"))
print(sampleTable)
fds = FraserDataSet(colData=sampleTable, workingDir=".", bamParam=ScanBamParam(mapqFilter=0), strandSpecific=0L)

getSplitReadCountsForAllSamples(fds)  # saves results to cache/
'""")
                    elif step == 2:
                        if sample_id in split_reads_jobs:
                            j.depends_on(split_reads_jobs[sample_id])
                        if j_extract_splice_junctions:
                            j.depends_on(j_extract_splice_junctions)

                        j.command(f"gsutil -m cp {output_file_path_splice_junctions_RDS} .")

                        # count non-split reads
                        j.command(f"""time xvfb-run Rscript -e '
library(FRASER)
library(data.table)

spliceJunctions = readRDS("{os.path.basename(output_file_path_splice_junctions_RDS)}")

sampleTable = data.table(sampleID=c("{sample_id}"), bamFile=c("{bam_path}"))
print(sampleTable)

fds = FraserDataSet(colData=sampleTable, workingDir=".", bamParam=ScanBamParam(mapqFilter=0), strandSpecific=0L)

getNonSplitReadCountsForAllSamples(fds, spliceJunctions)  # saves results to cache/
'""")
                    j.command(f"ls -lh .")
                    j.command(f"tar czf {os.path.basename(output_file_path)} cache")
                    j.command(f"gsutil -m cp {os.path.basename(output_file_path)} {output_file_path}")

                    j.command(f"echo Done: {output_file_path}")
                    j.command(f"date")

                    print("Output file path: ", output_file_path)

                    if step == 1:
                        split_reads_jobs[sample_id] = j
                    elif step == 2:
                        non_split_reads_jobs[sample_id] = j

                if len(split_reads_output_files) == 0:
                    break

                if step == 1 and not args.skip_step1:
                    if hl.hadoop_is_file(output_file_path_splice_junctions_RDS) and not args.force:
                        logger.info(f"{output_file_path_splice_junctions_RDS} file already exists. Skipping extractSpliceJunctions step...")
                        continue

                    j_extract_splice_junctions = batch_utils.init_job(batch, f"{sample_set_label}: Extract splice-junctions", disk_size=30, memory=60, image=DOCKER_IMAGE)
                    for j in split_reads_jobs.values():
                        j_extract_splice_junctions.depends_on(j)

                    extract_splice_junctions(
                        j_extract_splice_junctions,
                        split_reads_output_files,
                        args.cpu,
                        output_file_path_splice_junctions_RDS)

                elif step == 2 and not args.skip_step2:
                    if hl.hadoop_is_file(output_file_path_calculated_psi_values_tar_gz) and not args.force:
                        logger.info(f"{output_file_path_calculated_psi_values_tar_gz} file already exists. Skipping calculatePSIValues step...")
                        continue

                    num_cpu = 4 if args.local else 16
                    memory = 60
                    j_calculate_psi_values = batch_utils.init_job(batch, f"{sample_set_label}: Calculate PSI values", disk_size=50, cpu=num_cpu, memory=memory, image=DOCKER_IMAGE)
                    if j_extract_splice_junctions:
                        j_calculate_psi_values.depends_on(j_extract_splice_junctions)
                    for j in non_split_reads_jobs.values():
                        j_calculate_psi_values.depends_on(j)

                    calculate_psi_values(
                        j_calculate_psi_values,
                        sample_set_label,
                        split_reads_output_files,
                        non_split_reads_output_files,
                        output_file_path_splice_junctions_RDS,
                        args.metadata_tsv_path,
                        output_file_path_calculated_psi_values_tar_gz)

            # filter and annotate data
            if args.skip_step3:
                logger.info(f"Skipping filterAndAnnotate step...")
            elif hl.hadoop_is_file(output_file_path_filter_and_annotate_data_tar_gz) and not args.force:
                logger.info(f"{output_file_path_filter_and_annotate_data_tar_gz} file already exists. Skipping calculatedBestQ step...")
            else:
                num_cpu = 4 if args.local else 16
                memory = 3.75*num_cpu
                j_filter_and_annotate_data = batch_utils.init_job(batch, f"{sample_set_label}: Filter and Annotate", disk_size=50, cpu=num_cpu, memory=memory, image=DOCKER_IMAGE)

                if j_calculate_psi_values:
                    j_filter_and_annotate_data.depends_on(j_calculate_psi_values)

                filter_and_annotate_data(
                    j_filter_and_annotate_data,
                    sample_set_label,
                    output_file_path_calculated_psi_values_tar_gz,
                    output_file_path_filter_and_annotate_data_tar_gz)

            # compute Best Q
            if args.skip_step4:
                logger.info(f"Skipping calculatedBestQ step...")
            elif hl.hadoop_is_file(output_file_path_calculated_best_q_tar_gz) and not args.force:
                logger.info(f"{output_file_path_calculated_best_q_tar_gz} file already exists. Skipping calculatedBestQ step...")
            else:
                num_cpu = 4 if args.local else 16
                memory = 3.75*num_cpu

                j_calculate_best_q = batch_utils.init_job(batch, f"{sample_set_label}: Calculate Best Q", disk_size=50, cpu=num_cpu, memory=memory, image=DOCKER_IMAGE)

                if j_filter_and_annotate_data:
                    j_calculate_best_q.depends_on(j_filter_and_annotate_data)

                calculate_best_q(
                    j_calculate_best_q,
                    sample_set_label,
                    4,
                    output_file_path_filter_and_annotate_data_tar_gz,
                    output_file_path_calculated_best_q_tar_gz)

            # output_file_path_fraser_analysis_tar_gz
            if args.skip_step5:
                logger.info(f"Skipping FRASER analysis step...")
            elif hl.hadoop_is_file(output_file_path_fraser_analysis_results_only_tar_gz) and not args.force:
                logger.info(f"{output_file_path_fraser_analysis_results_only_tar_gz} file already exists. Skipping run_fraser_analysis step...")
            else:
                num_cpu = 4 if args.local else 16
                memory = 3.75*num_cpu

                j_fraser_analysis = batch_utils.init_job(batch, f"{sample_set_label}: Run Fraser Analysis", disk_size=50, cpu=num_cpu, memory=memory, image=DOCKER_IMAGE)
                if j_calculate_best_q:
                    j_fraser_analysis.depends_on(j_calculate_best_q)

                run_fraser_analysis(
                    j_fraser_analysis,
                    sample_set_label,
                    output_file_path_calculated_best_q_tar_gz,
                    output_file_path_fraser_analysis_tar_gz,
                    output_file_path_fraser_analysis_results_only_tar_gz)

            if hl.hadoop_is_file(output_file_path_fraser_volcano_plots_tar_gz) and not args.force:
                logger.info(f"{output_file_path_fraser_volcano_plots_tar_gz} file already exists. Skipping run_fraser_analysis step...")
            else:
                num_cpu = 4
                memory = 3.75*num_cpu

                j_get_plots = batch_utils.init_job(batch, f"{sample_set_label}: Get Plots", disk_size=50, cpu=num_cpu, memory=memory, image=DOCKER_IMAGE)
                if j_fraser_analysis:
                    j_get_plots.depends_on(j_fraser_analysis)

                run_plots(
                    j_get_plots,
                    sample_set_label,
                    output_file_path_fraser_analysis_tar_gz,
                    output_file_path_fraser_volcano_plots_tar_gz)


if __name__ == "__main__":
    main()
