import datetime
import hail as hl
import hashlib
import logging
import os
import pandas as pd

from batch import batch_utils
from sample_metadata.rnaseq_metadata_utils import get_rnaseq_downstream_analysis_metadata_df
from gagneurlab.gagneur_utils import GENCODE_TXDB, DOCKER_IMAGE, GCLOUD_PROJECT, GCLOUD_CREDENTIALS_LOCATION, GCLOUD_USER_ACCOUNT, OUTRIDER_COUNTS_TSV_GZ


logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s', level=logging.INFO)
logger = logging.getLogger(__name__)

hl.init(log="/dev/null")

POSSIBLE_CONFOUNDERS = """c("tissue", "sex", "stranded", "read_length", "batch")""" # "RIN"
PADJ_THRESHOLD = 0.05

NUM_CPU = 16
def main():
    downstream_analysis_df = get_rnaseq_downstream_analysis_metadata_df()

    p = batch_utils.init_arg_parser(default_cpu=NUM_CPU, gsa_key_file=os.path.expanduser("~/.config/gcloud/misc-270914-cb9992ec9b25.json"))
    p.add_argument("--counts-tsv-path", default=OUTRIDER_COUNTS_TSV_GZ, help="Counts .tsv")
    p.add_argument("--skip-step1", action="store_true", help="Skip initial steps including computing best Q")
    p.add_argument("--skip-step2", action="store_true", help="Skip OUTRIDER fit step")

    g = p.add_mutually_exclusive_group()
    g.add_argument("--with-gtex", help="Use GTEX controls.", action="store_true")
    g.add_argument("--only-gtex", help="Run on just the GTEX control samples to test FP rate.", action="store_true")

    p.add_argument("tissue", nargs="+", choices=set(downstream_analysis_df["tissue"]), help="Name of RNA-seq batch to process")
    args = p.parse_args()

    #local_all_counts_tsv_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.path.basename(OUTRIDER_COUNTS_TSV_GZ))
    #print(f"Copying {local_all_counts_tsv_path} to {OUTRIDER_COUNTS_TSV_GZ}")
    #hl.hadoop_copy(local_all_counts_tsv_path, OUTRIDER_COUNTS_TSV_GZ)

    # process samples
    batch_label = f"OUTRIDER"
    if args.with_gtex:
        batch_label += " (with GTEx)"
    batch_label += ": "
    batch_label += ','.join(args.tissue)
    with batch_utils.run_batch(args, batch_label) as batch:

        for tissue in args.tissue:
            df = downstream_analysis_df[downstream_analysis_df["tissue"] == tissue]
            if not args.with_gtex and not args.only_gtex:
                df = df[~df["sample_id"].str.startswith("GTEX")]
            elif args.only_gtex:
                df = df[df["sample_id"].str.startswith("GTEX")]

            sample_ids = list(sorted(set(df.sample_id)))
            byte_string = ", ".join(sample_ids).encode()
            h = hashlib.md5(byte_string).hexdigest().upper()
            sample_set_label = f"{tissue}__{len(sample_ids)}_samples_{h[:10]}"

            logger.info(f"Processing {sample_set_label}")

            c_vector_of_sample_names = 'c("' + '", "'.join(sample_ids) + '")'
            print(f"Samples: {c_vector_of_sample_names}")

            if args.with_gtex:
                sample_set_label += "_with_GTEX"
            elif args.only_gtex:
                sample_set_label += "_only_GTEX"
            else:
                sample_set_label += "_without_GTEX"

            output_base_dir = f"gs://tgg-rnaseq/gagneur/outrider/results/"
            if "sequencing_date" in set(df.columns):
                most_recent_sequencing_date = str(max(df.sequencing_date)).replace("-", "_")
                output_base_dir += f"{most_recent_sequencing_date}__"
            output_base_dir += f"{sample_set_label}"

            # upload metadata.tsv
            metadata_tsv_filename = f"sample_metadata_{sample_set_label}.tsv"
            local_metadata_tsv_path = f"/tmp/{metadata_tsv_filename}"
            metadata_tsv_path = os.path.join(output_base_dir, metadata_tsv_filename)
            df.to_csv(local_metadata_tsv_path, header=True, index=False, sep="\t")
            hl.hadoop_copy(local_metadata_tsv_path, metadata_tsv_path)

            j = batch_utils.init_job(batch, sample_set_label, image=DOCKER_IMAGE if not args.raw else None, cpu=NUM_CPU, memory="highmem", disk_size=50)
            batch_utils.switch_gcloud_auth_to_user_account(j, GCLOUD_CREDENTIALS_LOCATION, GCLOUD_USER_ACCOUNT, GCLOUD_PROJECT)
            # copy inputs
            j.command(f"""cd /io""")
            j.command(f"""gsutil -m cp {GENCODE_TXDB} .""")
            j.command(f"""gsutil -m cp {metadata_tsv_path} {args.counts_tsv_path} .""")
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
sample_info_path = "{os.path.basename(metadata_tsv_path)}"
sampleInfo = fread(sample_info_path)
sampleInfo$read_length = as.character(sampleInfo$read_length)

sampleSetLabel = "{sample_set_label}"
sampleSubset = {c_vector_of_sample_names}
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

#print("=============================================================")
#print("sampleInfo")
#print(sampleInfo, topn=100, nrows=1000)
#print("=============================================================")

cnts = cnts[,sampleInfo$sampleID]
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
    ods = findEncodingDim(ods, BPPARAM=MulticoreParam({NUM_CPU}, progressbar=TRUE))
    g = plotEncDimSearch(ods)
    ggsave(file=paste(sampleSetLabel, "__plotEncDimSearch", ".png", sep=""), g, type="cairo")
    optimal_q = metadata(ods)$opt
}} else {{
    optimal_q = length(sampleSubset)
    metadata(ods)$opt = optimal_q 
}}

saveRDS(ods, "{os.path.basename(step1_output_RDS_file)}")
'""")
                j.command(f"gsutil -m cp  *.tsv *.pdf *.png {output_base_dir}/")
                j.command(f"""gsutil -m cp "{os.path.basename(step1_output_RDS_file)}" {output_base_dir}/""")

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

ods = OUTRIDER(ods, verbose=TRUE, iterations=15, q=q, BPPARAM=MulticoreParam({NUM_CPU}, progressbar=TRUE))
rownames(ods) <- gsub("\\\\.[0-9]*(_[0-9]*)?.*$", "", rownames(ods))

saveRDS(ods, "{os.path.basename(step2_output_RDS_file)}")

plotCountCorHeatmap(ods, colGroups=possibleConfounders, normalized=TRUE, device="pdf", type="cairo", nRowCluster=1, nColCluster=1, main=paste("Count correlation heatmap q=", q, sep=""), filename=paste(sampleSetLabel, "__plotCountCorHeatmap_after_correction.pdf", sep=""))

plotCountGeneSampleHeatmap(ods, colGroups=possibleConfounders, normalized=TRUE, device="pdf", type="cairo", main=paste("Count Gene vs Sample Heatmap q=", q, sep=""), device="pdf", type="cairo", filename=paste(sampleSetLabel, "__plotCountGeneSampleHeatmap_after_correction.pdf", sep=""))

g = plotAberrantPerSample(ods, padjCutoff={PADJ_THRESHOLD})
ggsave(file=paste(sampleSetLabel, "__aberrantPerSample_padj_{PADJ_THRESHOLD}.png", sep=""), g, type="cairo")
'""")
                j.command(f"gsutil -m cp  *.tsv *.pdf *.png {output_base_dir}/")
                j.command(f"""gsutil -m cp "{os.path.basename(step2_output_RDS_file)}" {output_base_dir}/""")

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

resultTableColumns = c("sampleID", "geneID", "symbol", "biotype", "pValue", "padjust", "zScore", "rawcounts", "normcounts", "meanCorrected", "description", "chr", "start", "end", "strand")
q = metadata(ods)$opt

res = results(ods, all=TRUE)
res = merge(res, grch38, by.x="geneID", by.y="ensgene", sort=FALSE, all.x=TRUE)
res = res[, ..resultTableColumns]
setorderv(res, c("sampleID", "padjust"))
res[, "q"] = q
write.table(res, file=paste(sampleSetLabel, "__ods__", "q", q, "_all_results.tsv", sep=""), quote=FALSE, sep="\\t", row.names=FALSE)

res = results(ods, padjCutoff={PADJ_THRESHOLD})
res = merge(res, grch38, by.x="geneID", by.y="ensgene", sort=FALSE, all.x=TRUE)[!duplicated(geneID),]
res = res[, ..resultTableColumns]
setorderv(res, c("sampleID", "padjust"))
res[, "q"] = q
write.table(res, file=paste(sampleSetLabel, "__ods__", "q", q, "_padj_{PADJ_THRESHOLD}_results.tsv", sep=""), quote=FALSE, sep="\\t", row.names=FALSE)

# annotate gene names 
geneIdMap <- merge(data.table(ensgene=rownames(ods)), grch38, sort=FALSE, all.x=TRUE)[!duplicated(ensgene),]
geneLabels <- geneIdMap[, ifelse(is.na(symbol) | symbol == "" | duplicated(symbol), ensgene, symbol)]

sample_ids = unique(res$sampleID)
sample_ids = sample_ids[order(sample_ids)]
for(sample_id in sample_ids) {{
    print(paste("Plotting volcano plots for ", sample_id))
    volcanoPlotLabels = ifelse(names(ods) %in% res[res$sampleID == sample_id]$geneID, geneLabels, "")
    p = plotVolcano(ods, sample_id, padjCutoff={PADJ_THRESHOLD}, basePlot=TRUE) + geom_label_repel(aes(label=volcanoPlotLabels), force=3, nudge_y = -1, box.padding = 0.35, point.padding = 0.5, segment.color = "grey50") + theme(plot.margin=unit(c(0.5, 0, 0, 0), "cm"))
    ggsave(file=paste(sampleSetLabel, "__volcano__padj_{PADJ_THRESHOLD}_", sample_id, ".png", sep=""), p, width=12, height=8, device="png", type="cairo") 
}}
'""")
                j.command(f"""tar czf "{os.path.basename(step3_output_tar_gz_file)}" *.png *.tsv""")
                j.command(f"""gsutil -m cp "{os.path.basename(step3_output_tar_gz_file)}" {output_base_dir}/""")


if __name__ == "__main__":
    main()


