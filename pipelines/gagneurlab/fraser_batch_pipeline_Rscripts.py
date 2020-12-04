import os

PADJ_THRESHOLD=0.1
DELTA_PSI_THRESHOLD=0.3

def get_EXTRACT_SPLICE_JUNCTIONS_Rscript(bam_header_path, num_cpu):
    return f"""
library(FRASER)
library(data.table)
library(stringr)
library(purrr)
library(BiocParallel)

file_paths = list.files(".", pattern = "fraser_count_split_reads_.*.tar.gz$")
print(file_paths)
parse_sample_id = function(x) {{ return( str_replace(x[[1]], "fraser_count_split_reads_", "")) }}
sample_ids = unlist(map(strsplit(file_paths, "[.]"), parse_sample_id))

sampleTable = data.table(sampleID=sample_ids, bamFile="{os.path.basename(bam_header_path)}")
print(sampleTable)

if({num_cpu} == 1) {{
    bpparam = SerialParam(log=TRUE, progressbar=FALSE)
}} else {{
    bpparam = MulticoreParam({num_cpu}, log=FALSE, progressbar=FALSE)
}}

fds = FraserDataSet(colData=sampleTable, workingDir=".", bamParam=ScanBamParam(mapqFilter=0), strandSpecific=0L)
splitCountsForAllSamples = getSplitReadCountsForAllSamples(fds, BPPARAM=bpparam)
splitCountRanges = rowRanges(splitCountsForAllSamples)
print(splitCountRanges)

saveRDS(splitCountRanges, "spliceJunctions.RDS")
"""


def get_CALCULATE_PSI_VALUES_Rscript(splice_junctions_RDS_path, metadata_tsv_path, bam_header_path, num_cpu):
    return f"""
library(FRASER)
library(data.table)
library(stringr)
library(purrr)
library(BiocParallel)

splitCountRanges = readRDS("{os.path.basename(splice_junctions_RDS_path)}")
print(splitCountRanges)


file_paths = list.files(".", pattern = "fraser_count_split_reads_.*.tar.gz$")
print(file_paths)
parse_sample_id = function(x) {{ return( str_replace(x[[1]], "fraser_count_split_reads_", "")) }}
sample_ids = unlist(map(strsplit(file_paths, "[.]"), parse_sample_id))

sampleTable = fread("{os.path.basename(metadata_tsv_path)}")
sampleTable$read_length = as.character(sampleTable$read_length)

sampleTable = sampleTable[sampleTable$sample_id %in% sample_ids]
if (nrow(sampleTable) != length(sample_ids)) {{
    print(paste("ERROR: nrow(sampleTable) != length(sample_ids):", nrow(sampleTable), length(sample_ids)))
    quit("yes")
}}

sampleTable$bamFile = "{os.path.basename(bam_header_path)}"
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
nonSplitCountsForAllSamples = splitCountsForAllSamples = NULL
gc()

fds = calculatePSIValues(fds, BPPARAM=bpparam)
saveFraserDataSet(fds)
"""

def get_FILTER_AND_ANNOTATE_DATA_Rscript(sample_set_label):
    return f"""
library(FRASER)
library(data.table)
library(stringr)
library(purrr)
library(BiocParallel)
library(annotables)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(ggrepel)
library(plotly)
library(RColorBrewer)
library(ggsci)
library(ggplot2)
library(gtable)

fds = loadFraserDataSet(".")
sampleSetLabel = "{sample_set_label}"

options(repr.plot.width=5, repr.plot.height=4)
fds = filterExpressionAndVariability(fds, minDeltaPsi={DELTA_PSI_THRESHOLD}, minExpressionInOneSample=2, filter=FALSE)
g = plotFilterExpression(fds, bins=100)
ggsave(file=paste(sampleSetLabel, "_plotFilterExpression.png", sep=""), g, device="png", type="cairo")

g = plotFilterVariability(fds) + theme(legend.position="none")
ggsave(file=paste(sampleSetLabel, "_plotFilterVariability.png", sep=""), g, device="png", type="cairo")
g = NULL
gc()

print(paste(length(fds), "splice junctions before filtering"))
fds = fds[mcols(fds, type="j")[,"passed"],]
print(paste(length(fds), "splice junctions after filtering"))

fds = annotateRanges(fds, GRCh=38)
saveFraserDataSet(fds)
"""

def get_CALCULATE_BEST_Q_Rscript(sample_set_label, num_cpu):
    return f"""
library(FRASER)
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

fds = loadFraserDataSet(".")
if({num_cpu}L == 1L) {{
    bpparam = SerialParam(log=TRUE, progressbar=FALSE)
}} else {{
    bpparam = MulticoreParam({num_cpu}, log=FALSE, threshold = "DEBUG", progressbar=FALSE)
}}

sampleSetLabel = "{sample_set_label}"
for(param_type in c("psi5", "psi3", "psiSite")) {{
    print("===============")
    fds = optimHyperParams(fds, param_type, plot=FALSE, implementation="PCA", BPPARAM=bpparam)
    g = plotEncDimSearch(fds, type=param_type, plotType="auc") 
    ggsave(file=paste(sampleSetLabel, "_plotEncDimSearch_", param_type,"_AUC.png", sep=""), g, device="png", type="cairo")
    g = plotEncDimSearch(fds, type=param_type, plotType="loss") 
    ggsave(file=paste(sampleSetLabel, "_plotEncDimSearch_", param_type,"_loss.png", sep=""), g, device="png", type="cairo")
    
    print(paste(i, ": ", bestQ(fds, type=i), sep=""))
    gc()
}}

print("===============")
for(param_type in c("psi5", "psi3", "psiSite")) {{
    print(paste(param_type, bestQ(fds, type=param_type), sep=" "))
}}

saveFraserDataSet(fds)
"""


def get_RUN_FRASER_ANALYSIS_Rscript(sample_set_label, num_cpu):
    return f"""
library(FRASER)
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

fds = loadFraserDataSet(".")
if({num_cpu}L == 1L) {{
    bpparam = SerialParam(log=TRUE, progressbar=FALSE)
}} else {{
    bpparam = MulticoreParam({num_cpu}, log=FALSE, threshold = "DEBUG", progressbar=TRUE)
}}

sampleSetLabel = "{sample_set_label}"

possibleConfounders = c("tissue", "sex", "stranded", "read_length", "batch") 
for(param_type in c("psi5", "psi3", "psiSite")) {{
    plotCountCorHeatmap(fds, type=param_type, logit=TRUE, annotation_col=possibleConfounders, plotType="sampleCorrelation", device="pdf", filename=paste(sampleSetLabel, "_plotCountCorHeatmap_before_correction_", param_type ,".pdf", sep=""))
}}
for(param_type in c("psi5", "psi3", "psiSite")) {{
    plotCountCorHeatmap(fds, type=param_type, logit=TRUE, annotation_col=possibleConfounders, plotType="junctionSample", device="pdf", filename=paste(sampleSetLabel, "_plotCountJunctionSampleHeatmap_before_correction_", param_type ,".pdf", sep=""))
}}

gc()

message("Running FRASER with q=", bestQ(fds, type="psi5"), ", ", bestQ(fds, type="psi3"), " ", bestQ(fds, type="psiSite"))
implementation="PCA"
iterations = 15
for(param_type in c("psi5", "psi3", "psiSite")) {{
    q = bestQ(fds, type=param_type)
    message(date(), ": Fit step for: ", param_type, ". q=", q)
    fds = fit(fds, implementation = implementation, q = q, iterations=iterations, type=param_type, BPPARAM=bpparam)
    message(date(), ": Compute p values for: ", param_type, ".")
    fds = calculatePvalues(fds, type=param_type)
    message(date(), ": Adjust p values for: ", param_type, ".")
    fds = calculatePadjValues(fds, type=param_type)
    message(date(), ": Compute Z scores for: ", param_type, ".")
    fds = calculateZscore(fds, type = param_type)
  
    plot_filename1 = paste(sampleSetLabel, "_using", implementation, "_plotCountCorHeatmap_after_correction_", param_type, "__q", q, ".pdf", sep="")
    plot_filename2 = paste(sampleSetLabel, "_using", implementation, "_plotCountJunctionSampleHeatmap_after_correction_", param_type, "__q", q ,".pdf", sep="")
    message("Creating heatmap plot: ", plot_filename1)  
    plotCountCorHeatmap(fds, type=param_type, normalized=TRUE, logit=TRUE, annotation_col=possibleConfounders, plotType="sampleCorrelation", device="pdf", filename=plot_filename1)
    message("Creating heatmap plot: ", plot_filename2)  
    plotCountCorHeatmap(fds, type=param_type, normalized=TRUE, logit=TRUE, annotation_col=possibleConfounders, plotType="junctionSample", device="pdf", filename=plot_filename2)
    message("Done creating heatmap plots")  
    gc()
}}

qLabel = paste("_using", implementation, "_fds__psi5_q", bestQ(fds, type="psi5"), "__psi3_q", bestQ(fds, type="psi3"), "__psiSite_q", bestQ(fds, type="psiSite"), sep="")

g = plotAberrantPerSample(fds)
ggsave(file=paste(sampleSetLabel, "_plotAberrantPerSample.png", sep=""), g, device="png", type="cairo")

#res = results(fds, padjCutoff=1, zScoreCutoff=NA, deltaPsiCutoff=NA)
#saveRDS(res, paste(sampleSetLabel, qLabel, "_all_results.RDS", sep=""))
#message("Done saving results RDS") 
#gc()

#filename=paste(sampleSetLabel, qLabel, "_all_results.tsv.gz", sep="")
#write.table(as.data.table(res), file=filename, quote=FALSE, sep="\\t", row.names=FALSE)
#gc()

res = results(fds, padjCutoff={PADJ_THRESHOLD}, zScoreCutoff=NA, deltaPsiCutoff=NA)
saveRDS(res, paste(sampleSetLabel, qLabel, "_padj_{PADJ_THRESHOLD}_results.RDS", sep=""))
message("Done saving results RDS")  
gc()

message(length(res), " junctions in results with p < {PADJ_THRESHOLD}")
filename=paste(sampleSetLabel, qLabel, "_padj_{PADJ_THRESHOLD}_results.tsv.gz", sep="")
write.table(as.data.table(res), file=filename, quote=FALSE, sep="\\t", row.names=FALSE)
message("Done saving ", filename)
gc()

saveFraserDataSet(fds)

message("Done saving fds dataset")
"""


def get_PLOT_RESULTS(sample_set_label):
    return f"""
library(FRASER)
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

print("#### 1")
fds = loadFraserDataSet(".")
print("#### 2")
res = results(fds, padjCutoff={PADJ_THRESHOLD}, zScoreCutoff=NA, deltaPsiCutoff=NA)
print("#### 3")
sample_ids = res$sampleID[!duplicated(res$sampleID)]
print("#### 4")
sample_ids = sample_ids[order(sample_ids)]
print("#### 5")
for(sample_id in sample_ids) {{    
    print(paste("Plotting volcano plots for ", sample_id))
    for(param_type in c("psi5", "psi3", "psiSite")) {{
        print(paste("Plotting volcano plot", param_type, " for ", sample_id))
        volcanoPlotLabels = ifelse(names(fds) %in% res[sampleID == sample_id]$geneID, names(fds), "")
        p = plotVolcano(fds, type="psi5", sample_id, basePlot=TRUE) +
          geom_label_repel(aes(label=volcanoPlotLabels), force=3, nudge_y = -1, box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50') +
          labs(title=sample_id, x = "", y = "") +
          theme(plot.margin=unit(c(0.5, 0, 0, 0), "cm")) 
        
        ggsave(file=paste(sampleSetLabel, "_volcano_", param_type, "_", sample_id, ".png", sep=""), p, width=12, height=8, dpi=150)
    }}
}}

#plotQQ(fds, res[1, geneID], main="Q-Q plot for gene: CAPN3 @ q=20")
"""
