
# fraser
library(FRASER)
library(data.table)
library(ggplot2)
library(ggpubr)

#colab notebook: https://colab.research.google.com/drive/1OKT32eNIq7Cz839jjqz-GJlvoToPYbib

library(stringr)
library(dplyr)
library(purrr)
library(ggplot2)
library(data.table)
library(plotly)
library(stringr)
library(RColorBrewer)

library(ggsci)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)

txdb = TxDb.Hsapiens.UCSC.hg38.knownGene
orgDb = org.Hs.eg.db

# based on docs @ https://bioconductor.org/packages/devel/bioc/vignettes/FRASER/inst/doc/FRASER.pdf

setwd("~/project__rnaseq/data/samples/expression/bams")
bam_file_paths = list.files('.', pattern = "*.subset.bam$")
print(bam_file_paths)
sample_ids = unlist(map(strsplit(bam_file_paths, '[.]'), function(x) { return( str_replace(x[[1]], '.Aligned.sortedByCoord.out.subset.bam', '')) }))
sampleTable = data.table(sampleID=sample_ids, bamFile=bam_file_paths)

setwd("~/project__rnaseq/data/samples/fraser_output/analysis-node/gagneur/muscle/")
file_paths = list.files('cache/nonSplicedCounts/Data_Analysis/', pattern = "nonSplicedCounts-.*.h5$")
print(file_paths)
sample_ids = unlist(map(strsplit(file_paths, '[.]'), function(x) { return( str_replace(str_replace(x[[1]], 'nonSplicedCounts-', ''), '.h5$', '')) }))
sampleTable = data.table(sampleID=sample_ids, bamFile="~/project__rnaseq/data/samples/expression/bams/bam_header.bam") #ample_ids)

# https://github.com/c-mertes/FRASER/blob/master/R/FraserDataSet-class.R   - strandSpecific: 0 = no, 1 = yes, 2 = reverse
# strandSpecific = 1, 2 doesn't work for some


print(sampleTable)

fds = FraserDataSet(colData=sampleTable, workingDir=".", bamParam=ScanBamParam(mapqFilter=0), strandSpecific=0L)
splitCountsForAllSamples = getSplitReadCountsForAllSamples(fds)
splitCountRanges = rowRanges(splitCountsForAllSamples)
print(splitCountRanges)
splitCountRanges = readRDS("./spliceJunctions_188_samples_6A483623C9.RDS")
print(splitCountRanges)
nonSplitCountsForSample = getNonSplitReadCountsForAllSamples(fds, splitCountRanges)
#ds1 = addCountsToFraserDataSet(fds, splitCountsForSample, nonSplitCountsForSample)

fds = countRNAData(fds, NcpuPerSample=1)
fds = calculatePSIValues(fds)


# MUN_FAM1_CTRL_02
# MUN_FAM5_SIBLINGMDC1A_01_R1


fds = annotateRanges(fds, GRCh=38)

fds = readRDS('fdsWithPSIValues.RDS')

assayNames(fds)
# "rawCountsJ"            
# "psi5"
# "psi3"                   
# "psiSite"                
# "delta_psi5"             
# "delta_psi3"             
# "delta_psiSite"         

"
topN = 50000, topJ = 5000, minMedian = 1, main = NULL, normalized = FALSE, 
  show_rownames = FALSE, show_colnames = FALSE, minDeltaPsi = 0.1, 
  annotation_col = NA, annotation_row = NA, border_color = NA, 
  nClust = 5, plotType = c('sampleCorrelation', 'junctionSample'), 
  sampleClustering = NULL, plotMeanPsi = TRUE, plotCov = TRUE, 
"
plotCountCorHeatmap(fds, type="psi5", logit=TRUE)





sample_ids = unlist(map(strsplit(bam_file_paths_sibling, '[.]'), function(x) { return(x[[1]]) }))
sampleTable = data.table(sampleID=sample_ids, bamFile=bam_file_paths_sibling)
print(sampleTable)
fds2 = FraserDataSet(colData=sampleTable, workingDir=".", bamParam=ScanBamParam(mapqFilter=0), strandSpecific=0L)   
#name(fds1) = "fds2"
fds2 = countRNAData(fds, NcpuPerSample=1)
fds2 = calculatePSIValues(fds2)


sample_ids = unlist(map(strsplit(bam_file_paths, '[.]'), function(x) { return(x[[1]]) }))
sampleTable = data.table(sampleID=sample_ids, bamFile=bam_file_paths)
print(sampleTable)
fds_both = FraserDataSet(colData=sampleTable, junctions=attr(fds1, spliceSites=attr(fds1, 'nonSplicedReads'), workingDir=".", bamParam=ScanBamParam(mapqFilter=0), strandSpecific=0L)   
fds_both = countRNAData(fds_both, NcpuPerSample=1)
fds_both = calculatePSIValues(fds_both)




sample_ids = unlist(map(strsplit(bam_file_paths, '[.]'), function(x) { return(x[[1]]) }))
sampleTable = data.table(sampleID=sample_ids, bamFile=bam_file_paths)
print(sampleTable)
fds = FraserDataSet(colData=sampleTable, workingDir=".", bamParam=ScanBamParam(mapqFilter=0), strandSpecific=0L)   
fds = countRNAData(fds, NcpuPerSample=1)
fds = calculatePSIValues(fds)


# https://github.com/c-mertes/FRASER/blob/master/R/FraserDataSet-class.R   - strandSpecific: 0 = no, 1 = yes, 2 = reverse
# strandSpecific = 1, 2 doesn't work for some
fds1 = FraserDataSet(colData=sampleTable, workingDir=".", bamParam=ScanBamParam(mapqFilter=0), strandSpecific=0L)   
fds1 = countRNAData(fds, NcpuPerSample=1)


fds = calculatePSIValues(fds)

fds = filterExpressionAndVariability(fds, minExpressionInOneSample=5, minDeltaPsi=0.0, filter=FALSE)
plotFilterExpression(fds, bins=100)
plotCountCorHeatmap(fds, type="psi5", logit=TRUE, normalized=FALSE)
plotCountCorHeatmap(fds, type="psi5", logit=TRUE, normalized=FALSE, plotType="junctionSample", topJ=100, minDeltaPsi = 0.01)
fds = FRASER(fds, q=3)  # remove noise and fit model
plotCountCorHeatmap(fds, type="psi5", logit=TRUE, normalized=TRUE)
#fds = annotateRangesWithTxDb(fds, txdb=txdb, orgDb=orgDb)

sessionInfo()
res = results(fds)

plotVolcano(fds, type="psi3", "MUN_FAM4_ATYPICALMDC1A_01_R1")  # psi3, psi5, or psiSite  deltaPsiCutoff = 0.3, padjCutoff = 0.1,

# plotExpression
# plotExpectedVsObservedPsi
# plotAberrantPerSample
# plotQQ


library(FRASER)
library(data.table)
library(stringr)
library(purrr)

setwd("~/project__rnaseq/data/samples/fraser_output")
file_paths = list.files('.', pattern = "fraser_count_rna.*.tar.gz$")
print(file_paths)
parse_sample_id = function(x) { return( str_replace(x[[1]], 'fraser_count_rna_', '')) }
sample_ids = unlist(map(strsplit(file_paths, '[.]'), parse_sample_id))

sampleTable = data.table(sampleID=sample_ids, bamFile=".")
print(sampleTable)

fds = FraserDataSet(colData=sampleTable, workingDir=".", bamParam=ScanBamParam(mapqFilter=0), strandSpecific=0L)
splitCountsForAllSamples = getSplitReadCountsForAllSamples(fds)
splitCountRanges = rowRanges(splitCountsForAllSamples)
saveRDS(splitCountRanges, "spliceJunctions.RDS")




