
# outrider
# library(OUTRIDER)
# library(annotables)
# library(data.table)
# library(ggplot2)
# library(ggpubr)

# devtools::install_github('bw2/FRASER')

# BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
# BiocManager::install("org.Hs.eg.db")

# fraser
library(FRASER)
library(data.table)
library(ggplot2)
library(ggpubr)

#colab notebook: https://colab.research.google.com/drive/1OKT32eNIq7Cz839jjqz-GJlvoToPYbib

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

bam_file_paths = list.files('.', pattern = "*.bam$")
sample_ids = unlist(map(strsplit(bam_file_paths, '[.]'), function(x) { return(x[[1]]) }))

sampleTable = data.table(sampleID=sample_ids, bamFile=bam_file_paths)

# https://github.com/c-mertes/FRASER/blob/master/R/FraserDataSet-class.R   - strandSpecific: 0 = no, 1 = yes, 2 = reverse
# strandSpecific = 1, 2 doesn't work for some
fds = FraserDataSet(colData=sampleTable, workingDir=".", bamParam=ScanBamParam(mapqFilter=0), strandSpecific=0L)   
fds = countRNAData(fds, NcpuPerSample=6)
fds = calculatePSIValues(fds)
fds = filterExpressionAndVariability(fds, minExpressionInOneSample=20, minDeltaPsi=0.0, filter=FALSE)
plotFilterExpression(fds, bins=100)
plotCountCorHeatmap(fds, type="psi5", logit=TRUE, normalized=FALSE)
plotCountCorHeatmap(fds, type="psi5", logit=TRUE, normalized=FALSE, plotType="junctionSample", topJ=100, minDeltaPsi = 0.01)
fds = FRASER(fds, q=3)  # remove noise and fit model
plotCountCorHeatmap(fds, type="psi5", logit=TRUE, normalized=TRUE)
fds = annotateRangesWithTxDb(fds, txdb=txdb, orgDb=orgDb)

res = results(fds)

plotVolcano(fds, type="psi3", "MUN_FAM2_TOTALMDC1A1_02")  # psi3, psi5, or psiSite  deltaPsiCutoff = 0.3, padjCutoff = 0.1,

# plotExpression
# plotExpectedVsObservedPsi
# plotAberrantPerSample
# plotQQ



