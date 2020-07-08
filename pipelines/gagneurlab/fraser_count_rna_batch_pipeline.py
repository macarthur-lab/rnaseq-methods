import hail as hl
import hashlib
import logging
import os
import pandas as pd
import sys

from batch import batch_utils
from gagneurlab.gagneur_utils import GAGNEUR_BATCHES, ALL_METADATA_TSV, BAM_HEADER_PATH, GENCODE_TXDB, DOCKER_IMAGE, GCLOUD_PROJECT, GCLOUD_CREDENTIALS_LOCATION, GCLOUD_USER_ACCOUNT

logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s', level=logging.INFO)
logger = logging.getLogger(__name__)


def main():
    p = batch_utils.init_arg_parser(default_cpu=4, gsa_key_file=os.path.expanduser("~/.config/gcloud/misc-270914-cb9992ec9b25.json"))
    p.add_argument("--with-gtex", help="Use GTEX controls.", action="store_true")
    p.add_argument("--skip-step1", action="store_true", help="Skip count-split-reads step")
    p.add_argument("-m1", "--memory-step1", type=float, help="Batch: (optional) memory in gigabytes (eg. 3.75)", default=3.75)
    p.add_argument("-m2", "--memory-step2", type=float, help="Batch: (optional) memory in gigabytes (eg. 3.75)", default=3.75)
    p.add_argument("--metadata-tsv-path", default=ALL_METADATA_TSV, help="Table with columns: sample_id, bam_path, bai_path, batch")
    p.add_argument("batch_name", nargs="+", choices=GAGNEUR_BATCHES.keys(), help="Name of RNA-seq batch to process")
    args = p.parse_args()

    hl.init(log="/dev/null", quiet=True)

    with hl.hadoop_open(args.metadata_tsv_path) as f:
        samples_df_unmodified = pd.read_table(f).set_index("sample_id", drop=False)

    with batch_utils.run_batch(args) as batch:

        for batch_name in args.batch_name:
            samples_df = samples_df_unmodified
            batch_dict = GAGNEUR_BATCHES[batch_name]
            batch_tissue = batch_dict['tissue']
            batch_sex = batch_dict['sex']

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
            #j_calculate_psi_values = None

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
                    #input_bam = "gs://macarthurlab-rnaseq/temp/MUN_FAM5_SIBLINGMDC1A_01_R1.Aligned.sortedByCoord.out.subset.bam"
                    #input_bai = "gs://macarthurlab-rnaseq/temp/MUN_FAM5_SIBLINGMDC1A_01_R1.Aligned.sortedByCoord.out.subset.bam.bai"

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

                    if step == 1 and args.skip_step1:
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

                if step == 1:
                    output_file_path_splice_junctions_RDS = os.path.join(output_dir_for_batch_specific_data, f"spliceJunctions_{sample_set_label}.RDS")
                    if hl.hadoop_is_file(output_file_path_splice_junctions_RDS) and not args.force:
                        logger.info(f"{output_file_path_splice_junctions_RDS} file already exists. Skipping extractSpliceJunctions step...")
                        continue

                    j_extract_splice_junctions = batch_utils.init_job(batch, f"{sample_set_label}: Extract splice-junctions", disk_size=30, memory=60, image=DOCKER_IMAGE)
                    batch_utils.switch_gcloud_auth_to_user_account(j_extract_splice_junctions, GCLOUD_CREDENTIALS_LOCATION, GCLOUD_USER_ACCOUNT, GCLOUD_PROJECT)

                    for j in split_reads_jobs.values():
                        j_extract_splice_junctions.depends_on(j)

                    for split_reads_output_files_batch in [split_reads_output_files[i:i+10] for i in range(0, len(split_reads_output_files), 10)]:
                        j_extract_splice_junctions.command(f"gsutil -m cp {' '.join(split_reads_output_files_batch)} .")
                    j_extract_splice_junctions.command(f"gsutil -m cp {BAM_HEADER_PATH} .")
                    j_extract_splice_junctions.command(f"for i in fraser_count_split_reads*.tar.gz; do tar xzf $i; done")
                    j_extract_splice_junctions.command(f"pwd && ls -lh && date && echo ------- && find cache -name '*.*'")
                    j_extract_splice_junctions.command(f"""time xvfb-run Rscript -e '
library(FRASER)
library(data.table)
library(stringr)
library(purrr)
library(BiocParallel)

file_paths = list.files(".", pattern = "fraser_count_split_reads_.*.tar.gz$")
print(file_paths)
parse_sample_id = function(x) {{ return( str_replace(x[[1]], "fraser_count_split_reads_", "")) }}
sample_ids = unlist(map(strsplit(file_paths, "[.]"), parse_sample_id))

sampleTable = data.table(sampleID=sample_ids, bamFile="{os.path.basename(BAM_HEADER_PATH)}")
print(sampleTable)

if({args.cpu} == 1) {{
    bpparam = SerialParam(log=TRUE, progressbar=FALSE)
}} else {{
    bpparam = MulticoreParam({args.cpu}, log=FALSE, progressbar=FALSE)
}}

fds = FraserDataSet(colData=sampleTable, workingDir=".", bamParam=ScanBamParam(mapqFilter=0), strandSpecific=0L)
splitCountsForAllSamples = getSplitReadCountsForAllSamples(fds, BPPARAM=bpparam)
splitCountRanges = rowRanges(splitCountsForAllSamples)
print(splitCountRanges)

saveRDS(splitCountRanges, "spliceJunctions.RDS")
'""")
                    j_extract_splice_junctions.command(f"ls -lh .")
                    j_extract_splice_junctions.command(f"echo ===============; echo cache dir; echo ===============; find cache")
                    j_extract_splice_junctions.command(f"echo ===============; echo savedObects dir; echo ===============; find savedObjects")

                    j_extract_splice_junctions.command(f"gsutil -m cp spliceJunctions.RDS {output_file_path_splice_junctions_RDS}")

                    print("Output file path: ", output_file_path_splice_junctions_RDS)
                elif step == 2:
                    output_file_path = os.path.join(output_dir_for_batch_specific_data, f"calculatedPSIValues_{sample_set_label}.tar.gz")
                    if hl.hadoop_is_file(output_file_path) and not args.force:
                        logger.info(f"{output_file_path} file already exists. Skipping calculatePSIValues step...")
                        #continue
                    num_cpu = 4 if args.local else 16
                    memory = 60
                    j_calculate_psi_values = batch_utils.init_job(batch, f"{sample_set_label}: Calculate PSI values", disk_size=50, cpu=num_cpu, memory=memory, image=DOCKER_IMAGE)
                    batch_utils.switch_gcloud_auth_to_user_account(j_calculate_psi_values, GCLOUD_CREDENTIALS_LOCATION, GCLOUD_USER_ACCOUNT, GCLOUD_PROJECT)

                    if j_extract_splice_junctions:
                        j_calculate_psi_values.depends_on(j_extract_splice_junctions)
                    for j in non_split_reads_jobs.values():
                        j_calculate_psi_values.depends_on(j)

                    j_calculate_psi_values.command(f"mkdir -p /tmp/fraser/{sample_set_label}")  # work-around for https://github.com/c-mertes/FRASER/issues/11
                    j_calculate_psi_values.command(f"cd /tmp/fraser/{sample_set_label}")
                    for split_reads_output_files_batch in [split_reads_output_files[i:i+10] for i in range(0, len(split_reads_output_files), 10)]:
                        j_calculate_psi_values.command(f"gsutil -m cp {' '.join(split_reads_output_files_batch)} .")
                    for non_split_reads_output_files_batch in [non_split_reads_output_files[i:i+10] for i in range(0, len(non_split_reads_output_files), 10)]:
                        j_calculate_psi_values.command(f"gsutil -m cp {' '.join(non_split_reads_output_files_batch)} .")
                    j_calculate_psi_values.command(f"gsutil -m cp {output_file_path_splice_junctions_RDS} .")
                    j_calculate_psi_values.command(f"gsutil -m cp {args.metadata_tsv_path} .")

                    j_calculate_psi_values.command(f"gsutil -m cp {BAM_HEADER_PATH} .")

                    j_calculate_psi_values.command(f"for i in fraser_count_split_reads*.tar.gz; do tar xzf $i; done")
                    j_calculate_psi_values.command(f"for i in fraser_count_non_split_reads*.tar.gz; do tar xzf $i; done")
                    j_calculate_psi_values.command(f"rm cache/nonSplicedCounts/Data_Analysis/spliceSiteCoordinates.RDS")
                    j_calculate_psi_values.command(f"pwd && ls -lh && date && echo ------- && find cache -name '*.*'")
                    j_calculate_psi_values.command(f"""time xvfb-run Rscript -e '
library(FRASER)
library(data.table)
library(stringr)
library(purrr)
library(BiocParallel)

splitCountRanges = readRDS("{os.path.basename(output_file_path_splice_junctions_RDS)}")
print(splitCountRanges)


file_paths = list.files(".", pattern = "fraser_count_split_reads_.*.tar.gz$")
print(file_paths)
parse_sample_id = function(x) {{ return( str_replace(x[[1]], "fraser_count_split_reads_", "")) }}
sample_ids = unlist(map(strsplit(file_paths, "[.]"), parse_sample_id))

sampleTable = fread("{os.path.basename(args.metadata_tsv_path)}")
sampleTable$read_length = as.character(sampleTable$read_length)

sampleTable = sampleTable[sampleTable$sample_id %in% sample_ids]
if (nrow(sampleTable) != length(sample_ids)) {{
    print(paste("ERROR: nrow(sampleTable) != length(sample_ids):", nrow(sampleTable), length(sample_ids)))
    quit("yes")
}}

sampleTable$bamFile = "{os.path.basename(BAM_HEADER_PATH)}"
setnames(sampleTable, "sample_id", "sampleID")

fds = FraserDataSet(colData=sampleTable, workingDir=".", bamParam=ScanBamParam(mapqFilter=0), strandSpecific=0L)
if({num_cpu}L == 1L) {{
    bpparam = SerialParam(log=TRUE, progressbar=FALSE)
}} else {{
    bpparam = MulticoreParam({num_cpu}, log=FALSE, threshold = "DEBUG", progressbar=FALSE)
}}

splitCountsForAllSamples = getSplitReadCountsForAllSamples(fds, BPPARAM=bpparam)
nonSplitCountsForAllSamples = getNonSplitReadCountsForAllSamples(fds, splitCountRanges, BPPARAM=bpparam)
fds = addCountsToFraserDataSet(fds, splitCountsForAllSamples, nonSplitCountsForAllSamples)
fds = calculatePSIValues(fds, BPPARAM=bpparam)
fds = filterExpressionAndVariability(fds, minDeltaPsi=0.0, filter=FALSE)

fds = optimHyperParams(fds, "psi5", plot=FALSE, BPPARAM=bpparam)
print(paste("psi5:", bestQ(fds, type="psi5"), sep=" "))
fds = optimHyperParams(fds, "psi3", plot=FALSE, BPPARAM=bpparam)
print(paste("psi3:", bestQ(fds, type="psi3"), sep=" "))
fds = optimHyperParams(fds, "psiSite", plot=FALSE, BPPARAM=bpparam)
print(paste("psiSite:", bestQ(fds, type="psiSite"), sep=" "))

print("===============")
print(paste("psi5:", bestQ(fds, type="psi5"), sep=" "))
print(paste("psi3:", bestQ(fds, type="psi3"), sep=" "))
print(paste("psiSite:", bestQ(fds, type="psiSite"), sep=" "))

saveFraserDataSet(fds)
'""")
                    # #fds = annotateRanges(fds, GRCh=38)
                    j_calculate_psi_values.command(f"pwd")
                    j_calculate_psi_values.command(f"ls -lh .")
                    j_calculate_psi_values.command(f"echo ===============; echo cache dir; echo ===============; find cache")
                    j_calculate_psi_values.command(f"echo ===============; echo savedObects dir; echo ===============; find savedObjects")
                    j_calculate_psi_values.command(f"rm *.tar.gz *.bam")
                    j_calculate_psi_values.command(f"cd ..")
                    j_calculate_psi_values.command(f"tar czf {os.path.basename(output_file_path)} {sample_set_label}")
                    j_calculate_psi_values.command(f"gsutil -m cp {os.path.basename(output_file_path)} {output_file_path}")
                    print("Output file path: ", output_file_path)


if __name__ == "__main__":
    main()
