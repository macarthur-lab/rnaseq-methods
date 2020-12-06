import datetime
import hail as hl
import hashlib
import logging
import os
import pandas as pd

from batch import batch_utils
from sample_metadata.rnaseq_metadata_utils import ANALYSIS_BATCHES
from gagneurlab.gagneur_utils import ALL_METADATA_TSV, GENCODE_TXDB, DOCKER_IMAGE, GCLOUD_PROJECT, GCLOUD_CREDENTIALS_LOCATION, GCLOUD_USER_ACCOUNT

logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s', level=logging.INFO)
logger = logging.getLogger(__name__)

OUTRIDER_COUNTS_TSV_GZ = "gs://macarthurlab-rnaseq/gagneur/outrider/OUTRIDER_input_table_RDG_and_GTEX_counts_for_all_tissues.tsv.gz"

POSSIBLE_CONFOUNDERS = """c("tissue", "sex", "stranded", "read_length", "batch")""" # "RIN"
PADJ_THRESHOLD = 0.05

def main():
    p = batch_utils.init_arg_parser(default_cpu=16, gsa_key_file=os.path.expanduser("~/.config/gcloud/misc-270914-cb9992ec9b25.json"))
    p.add_argument("--metadata-tsv-path", default=ALL_METADATA_TSV, help="Table with columns: sample_id, bam_path, bai_path, batch")
    p.add_argument("--counts-tsv-path", default=OUTRIDER_COUNTS_TSV_GZ, help="Counts .tsv")
    p.add_argument("--skip-step1", action="store_true", help="Skip initial steps including computing best Q")
    p.add_argument("--skip-step2", action="store_true", help="Skip OUTRIDER fit step")

    g = p.add_mutually_exclusive_group()
    g.add_argument("--with-gtex", help="Use GTEX controls.", action="store_true")
    g.add_argument("--only-gtex", help="Run on just the GTEX control samples to test FP rate.", action="store_true")

    p.add_argument("batch_name", nargs="+", choices=ANALYSIS_BATCHES.keys(), help="Name of RNA-seq batch to process")
    args = p.parse_args()

    if not args.force:
        hl.init(log="/dev/null", quiet=True)

    local_metadata_tsv_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.path.basename(ALL_METADATA_TSV))
    metadata_tsv_df = pd.read_table(local_metadata_tsv_path)

    print(f"Copying {local_metadata_tsv_path} to {ALL_METADATA_TSV}")
    hl.hadoop_copy(local_metadata_tsv_path, ALL_METADATA_TSV)

    #local_all_counts_tsv_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.path.basename(OUTRIDER_COUNTS_TSV_GZ))
    #print(f"Copying {local_all_counts_tsv_path} to {OUTRIDER_COUNTS_TSV_GZ}")
    #hl.hadoop_copy(local_all_counts_tsv_path, OUTRIDER_COUNTS_TSV_GZ)

    # process samples
    batch_label = f"OUTRIDER"
    if args.with_gtex:
        batch_label += " (with GTEx)"
    batch_label += ": "
    batch_label += ','.join(args.batch_name)
    with batch_utils.run_batch(args, batch_label) as batch:

        for batch_name in args.batch_name:
            batch_dict = ANALYSIS_BATCHES[batch_name]
            batch_tissue = batch_dict['tissue']
            batch_sex = batch_dict['sex']

            batch_df = metadata_tsv_df[metadata_tsv_df.sample_id.isin(batch_dict['samples'])]
            byte_string = ", ".join(sorted(batch_df.sample_id)).encode()
            h = hashlib.md5(byte_string).hexdigest().upper()
            sample_set_label = f"{batch_name}__{len(batch_df.sample_id)}_samples_{h[:10]}"

            num_sample_missing_from_metadata_tsv = len(batch_dict['samples']) - len(batch_df)
            if num_sample_missing_from_metadata_tsv > 0:
                raise ValueError(f"{batch_tissue} metadata_tsv_df is missing {num_sample_missing_from_metadata_tsv} sample ids. Download count files and run export_gagneur_metadata_table.py")
            logger.info(f"Processing {sample_set_label}")

            c_vector_of_sample_names = 'c("' + '", "'.join(batch_dict['samples']) + '")'
            if args.with_gtex:
                batch_include_GTEX_samples = "TRUE"
                sample_set_label += "_with_GTEX"
            elif args.only_gtex:
                c_vector_of_sample_names = "c()"
                batch_include_GTEX_samples = "TRUE"
                sample_set_label += "_only_GTEX"
            else:
                batch_include_GTEX_samples = "FALSE"
                sample_set_label += "_without_GTEX"

            output_base_dir = f"gs://macarthurlab-rnaseq/gagneur/outrider/results/{sample_set_label}/"

            j = batch_utils.init_job(batch, sample_set_label, DOCKER_IMAGE if not args.raw else None, args.cpu, args.cpu*3.75)
            batch_utils.switch_gcloud_auth_to_user_account(j, GCLOUD_CREDENTIALS_LOCATION, GCLOUD_USER_ACCOUNT, GCLOUD_PROJECT)
            # copy inputs
            j.command(f"""gsutil -m cp {GENCODE_TXDB} .""")
            j.command(f"""gsutil -m cp {args.metadata_tsv_path} {args.counts_tsv_path} .""")
            step1_output_RDS_file = os.path.join(output_base_dir, f"{sample_set_label}__ods.RDS")
            step2_output_RDS_file = os.path.join(output_base_dir, f"{sample_set_label}__ods_after_finding_bestQ.RDS")
            step3_output_tar_gz_file = os.path.join(output_base_dir, f"{sample_set_label}__volcano_plots__padj_{PADJ_THRESHOLD}.tar.gz")

            if args.skip_step1:
                logger.info(f"Skipping step 1...")
            elif not args.force and hl.hadoop_is_file(step1_output_RDS_file):
                logger.info(f"Step 1: output file exists: {step1_output_RDS_file} . Skipping step 1...")
            else:
                logger.info(f"Step 1: output file: {step1_output_RDS_file}")
                j.command(f"""time xvfb-run Rscript -e '

# outrider 
library(OUTRIDER)
library(annotables)
library(data.table)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(purrr)
library(ggrepel)
library(plotly)
library(stringr)
library(RColorBrewer)
library(ggsci)
library(ggplot2)
library(gtable)
library(grid)
library(gridExtra)

possibleConfounders = {POSSIBLE_CONFOUNDERS} 

# input tables generated by ~/project__rnaseq/code/rnaseq_methods/pipelines/gagneurlab/metadata/export_gagneur_metadata_table.py
# batches generated by ~/project__rnaseq/code/rnaseq_methods/pipelines/gagneurlab/metadata/metadata_notebook.py
sample_info_path = "{os.path.basename(args.metadata_tsv_path)}"
sampleInfo = fread(sample_info_path)
sampleInfo$read_length = as.character(sampleInfo$read_length)

GTEX_sampleIds = c()
if ({batch_include_GTEX_samples}) {{
    if (("{batch_sex}" == "M") || ("{batch_sex}" == "F")) {{
        GTEX_sampleIds = sampleInfo[(sampleInfo$sex == "{batch_sex}") & (sampleInfo$tissue == "{batch_tissue}") & grepl("GTEX", sampleInfo$sample_id)]$sample_id
    }} else {{
        GTEX_sampleIds = sampleInfo[(sampleInfo$tissue == "{batch_tissue}") & grepl("GTEX", sampleInfo$sample_id)]$sample_id    
    }}
}}


sampleSetLabel = "{sample_set_label}"
sampleSubset = {c_vector_of_sample_names}
sampleSubset = c(sampleSubset, GTEX_sampleIds)
print("sampleSubset: ")
print(sampleSubset)

sampleInfo = sampleInfo[sampleInfo$sample_id %in% sampleSubset]
print("sampleInfo$sample_id: ")
print(sampleInfo$sample_id)

if (length(sampleInfo$sample_id) != length(sampleSubset)) {{
    print(paste("ERROR: length(sampleInfo$sample_id) != length(sampleSubset):", length(sampleInfo$sample_id), length(sampleSubset)))
    print(paste("sampleInfo is missing:", setdiff(sampleSubset, sampleInfo$sample_id)))
    quit("yes")
}}

geneReadCounts = fread("{os.path.basename(args.counts_tsv_path)}", select=c("gene_id", sampleSubset))
geneReadCounts = geneReadCounts[!grep("ERCC", geneReadCounts$geneId),]


geneIds = geneReadCounts$gene_id
colsMiusGeneId = colnames(geneReadCounts)[!colnames(geneReadCounts) %in% c("gene_id")]
geneReadCounts = geneReadCounts[,..colsMiusGeneId]
rownames(geneReadCounts) = geneIds

cnts = as.matrix(geneReadCounts)
rownames(cnts) = geneIds
ncol(cnts)
nrow(cnts)
if (ncol(cnts) != length(sampleSubset)) {{
    print(paste("ERROR: ncol(cnts) != length(sampleSubset):", ncol(cnts), length(sampleSubset)))
    quit("yes")
}}

sampleInfo[,sampleID:=sample_id]
ods <- OutriderDataSet(countData=cnts, colData=sampleInfo)

ods <- estimateSizeFactors(ods)
sortedSizeFactors = sort(sizeFactors(ods))
g = ggplot(data=NULL, aes(y=sortedSizeFactors, x=1:ncol(ods))) + 
  geom_point(color="blue", size=1) + 
  labs(x="Sample rank", y="Size factors", title="Size factor distribution") + 
  geom_label_repel(aes(label=ifelse(sortedSizeFactors > 1.5, names(sortedSizeFactors), "")), 
                   nudge_x = -35, box.padding = 0.35, point.padding = 0.5, segment.color = "grey50") +
  geom_label_repel(aes(label=ifelse(sortedSizeFactors < 0.5, names(sortedSizeFactors), "")), 
                   nudge_x = 35, box.padding   = 0.35, point.padding = 0.5, segment.color = "grey50") +
  theme_bw()

ggsave(file=paste(sampleSetLabel, "__sizeFactors.png", sep=""), g, type="cairo")

print(sort(sizeFactors(ods))[1:5])

txdb <- loadDb("{os.path.basename(GENCODE_TXDB)}")
ods <- filterExpression(ods, gtfFile=txdb, filterGenes=FALSE)   #, fpkmCutoff=100)

g = plotFPKM(ods) + theme_bw() + theme(legend.position="bottom")
ggsave(file=paste(sampleSetLabel, "__plotFPKM.png", sep=""), g, device="png", type="cairo")

#plotExpressedGenes(ods)

print(paste(length(ods), "genes before filtering"))
ods <- ods[mcols(ods)$passedFilter,]
print(paste(length(ods), "genes after filtering"))

plotCountCorHeatmap(ods, colGroups=possibleConfounders, normalized=FALSE, device="pdf", type="cairo", nRowCluster=1, nColCluster=1, filename=paste(sampleSetLabel, "__plotCountCorHeatmap_before_correction.pdf", sep=""))
plotCountGeneSampleHeatmap(ods, colGroups=possibleConfounders, normalized=FALSE, device="pdf", type="cairo", filename=paste(sampleSetLabel, "__plotCountGeneSampleHeatmap_before_correction.pdf", sep=""))

if (length(sampleSubset) > 5) {{
    ods = findEncodingDim(ods, BPPARAM=MulticoreParam({args.cpu}, progressbar=TRUE))
    g = plotEncDimSearch(ods)
    ggsave(file=paste(sampleSetLabel, "__plotEncDimSearch", ".png", sep=""), g, type="cairo")
    optimal_q = metadata(ods)$opt
}} else {{
    optimal_q = length(sampleSubset)
    metadata(ods)$opt = optimal_q 
}}

saveRDS(ods, "{os.path.basename(step1_output_RDS_file)}")
'""")
                j.command(f"gsutil -m cp  *.pdf *.png {output_base_dir}")
                j.command(f"""gsutil -m cp "{os.path.basename(step1_output_RDS_file)}" {output_base_dir}""")

            if args.skip_step2:
                logger.info(f"Skipping step 2...")
            elif not args.force and hl.hadoop_is_file(step2_output_RDS_file):
                logger.info(f"Step 2: output file exists: {step2_output_RDS_file} . Skipping step 2...")
            else:
                logger.info(f"Step 2: output file: {step2_output_RDS_file}")

                j.command(f"gsutil -m cp -n {step1_output_RDS_file} .")

                j.command(f"""time xvfb-run Rscript -e '

library(OUTRIDER)
library(annotables)
library(data.table)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(purrr)
library(ggrepel)
library(plotly)
library(stringr)
library(RColorBrewer)
library(ggsci)
library(ggplot2)
library(gtable)
library(grid)
library(gridExtra)

sampleSetLabel = "{sample_set_label}"
possibleConfounders = {POSSIBLE_CONFOUNDERS} 

ods = readRDS("{os.path.basename(step1_output_RDS_file)}")

q = metadata(ods)$opt

ods = OUTRIDER(ods, verbose=TRUE, iterations=15, q=q, BPPARAM=MulticoreParam({args.cpu}, progressbar=TRUE))
saveRDS(ods, "{os.path.basename(step2_output_RDS_file)}")

plotCountCorHeatmap(ods, colGroups=possibleConfounders, normalized=TRUE, device="pdf", type="cairo", nRowCluster=1, nColCluster=1, main=paste("Count correlation heatmap q=", q, sep=""), filename=paste(sampleSetLabel, "__plotCountCorHeatmap_after_correction.pdf", sep=""))

plotCountGeneSampleHeatmap(ods, colGroups=possibleConfounders, normalized=TRUE, device="pdf", type="cairo", main=paste("Count Gene vs Sample Heatmap q=", q, sep=""), device="pdf", type="cairo", filename=paste(sampleSetLabel, "__plotCountGeneSampleHeatmap_after_correction.pdf", sep=""))

g = plotAberrantPerSample(ods, padjCutoff={PADJ_THRESHOLD})
ggsave(file=paste(sampleSetLabel, "__aberrantPerSample_padj_{PADJ_THRESHOLD}.png", sep=""), g, type="cairo")

# annotate gene names  
geneIDs <- gsub("\\\\.[0-9]*(_[0-9]*)?.*$", "", rownames(ods))
geneIdMap <- merge(data.table(ensgene=geneIDs), grch38, sort=FALSE, all.x=TRUE)[!duplicated(ensgene),]
  
# set new gene names only if hgnc symbol is present
if(!"ENSG" %in% colnames(mcols(ods))) {{
    mcols(ods)$ENSG <- geneIDs
    rownames(ods) <- geneIdMap[,ifelse(is.na(symbol) | symbol == "" | duplicated(symbol), geneIDs, symbol)]
}}

res = results(ods, padjCutoff=1)
res = res[,c("sampleID", "geneID", "pValue", "padjust", "zScore", "rawcounts")][order(padjust),]
res[, "q"] = q
write.table(res, file=paste(sampleSetLabel, "__ods__", "q", q, "_all_results.tsv.gz", sep=""), quote=FALSE, sep="\\t", row.names=FALSE)

res = results(ods, padjCutoff={PADJ_THRESHOLD})
res = res[,c("sampleID", "geneID", "pValue", "padjust", "zScore", "rawcounts")][order(padjust),]
res[, "q"] = q
write.table(res, file=paste(sampleSetLabel, "__ods__", "q", q, "_padj_{PADJ_THRESHOLD}_results.tsv.gz", sep=""), quote=FALSE, sep="\\t", row.names=FALSE)
'""")
                j.command(f"gsutil -m cp  *.tsv.gz *.pdf *.png {output_base_dir}")
                j.command(f"""gsutil -m cp "{os.path.basename(step2_output_RDS_file)}" {output_base_dir}""")

            if not args.force and hl.hadoop_is_file(step3_output_tar_gz_file):
                logger.info(f"Step 3: output file exists: {step3_output_tar_gz_file} . Skipping step 3...")
            else:
                logger.info(f"Step 3: output file: {step3_output_tar_gz_file}")

                j.command(f"gsutil -m cp -n {step2_output_RDS_file} .")
                j.command(f"""time xvfb-run Rscript -e '
library(OUTRIDER)
library(annotables)
library(data.table)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(purrr)
library(ggrepel)
library(plotly)
library(stringr)
library(RColorBrewer)
library(ggsci)
library(ggplot2)
library(gtable)
library(grid)
library(gridExtra)

sampleSetLabel = "{sample_set_label}"
ods = readRDS("{os.path.basename(step2_output_RDS_file)}")

res = results(ods)

sample_ids = unique(res$sampleID)
sample_ids = sample_ids[order(sample_ids)]
for(sample_id in sample_ids) {{
    print(paste("Plotting volcano plots for ", sample_id))
    volcanoPlotLabels = ifelse(names(ods) %in% res[res$sampleID == sample_id]$geneID, names(ods), "")
    p = plotVolcano(ods, sample_id, padjCutoff={PADJ_THRESHOLD}, basePlot=TRUE) + geom_label_repel(aes(label=volcanoPlotLabels), force=3, nudge_y = -1, box.padding = 0.35, point.padding = 0.5, segment.color = "grey50") + theme(plot.margin=unit(c(0.5, 0, 0, 0), "cm"))
    ggsave(file=paste(sampleSetLabel, "__volcano__padj_{PADJ_THRESHOLD}_", sample_id, ".png", sep=""), p, width=12, height=8, device="png", type="cairo") 
}}
'""")
                j.command(f"""tar czf "{os.path.basename(step3_output_tar_gz_file)}" *.png""")
                j.command(f"""gsutil -m cp "{os.path.basename(step3_output_tar_gz_file)}" {output_base_dir}""")


if __name__ == "__main__":
    main()


