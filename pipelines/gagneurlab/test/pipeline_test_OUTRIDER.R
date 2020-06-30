# colab notebook: https://colab.research.google.com/drive/1OKT32eNIq7Cz839jjqz-GJlvoToPYbib

# outrider 
library(OUTRIDER)
library(annotables)
library(data.table)
library(ggplot2)
library(ggpubr)

# fraser
#library(FraseR)
#library(data.table)
#library(ggplot2)
#library(ggpubr)



library(dplyr)
library(purrr)
library(ggrepel)
library(data.table)
library(plotly)
library(stringr)
library(RColorBrewer)
library(ggsci)
library(ggplot2)
library(gtable)
library(grid)
library(gridExtra)

possibleConfounders = c("tissue", "sex", "stranded", "read_length", 'batch') # "RIN"')

sander_muscle = c('HK069-0177_2', 'HK072-001_2', 'OUN_HK018_0047', 'OUN_HK047_0116', 'OUN_HK079_001', 'OUN_HK080_001', 'OUN_HK081_001', 'OUN_HK112_001', 'OUN_HK116_001', 'OUN_HK124_001', 'OUN_HK137_001')
#sander_muscle = c('OUN_HK079_001')
sander_muscle = c('HK072-001_2')


sander_fibs = c('HK006_0016', 'HK010_0026', 'HK011_0029', 'HK024_0065', 'HK044_0109', 'HK085-001', 'HK088-001_2', 'HK101-001', 'HK106-001')
#RGP_muscle = c('RGP_248_3', 'RGP_273_3', 'RGP_54_3_2', 'RGP_56_3_3', 'RGP_800_3_2', 'RGP_552_3')
#RGP_fibs = c('RGP_7_1_2',  'RGP_7_2_2', 'RGP_7_3_2', 'RGP_7_4_2', 'RGP_7_5_2')

#RGP_blood = c('RGP_94_3_2', 'RGP_94_4_2', 'RGP_554_3_2')


setwd("~/project__rnaseq/code/rnaseq_methods/pipelines/gagneurlab/metadata")
#sampleInfo = fread("metadata_table_for_whole_blood.tsv")
#geneReadCounts = fread("OUTRIDER_input_table_for_whole_blood.tsv")
#sampleInfo = fread("metadata_table_for_lymphocytes.tsv")
#geneReadCounts = fread("OUTRIDER_input_table_for_lymphocytes.tsv")
#sampleInfo = fread("metadata_table_for_fibroblasts.tsv")
#geneReadCounts = fread("OUTRIDER_input_table_for_fibroblasts.tsv")
sampleInfo = fread("metadata_table_for_muscle.tsv")
GTEX_muscle_sample_ids = sampleInfo$sample_id[grep("GTEX", sampleInfo$sample_id)]
#geneReadCounts = fread("OUTRIDER_input_table_for_muscle.tsv")
#sampleInfo = fread("metadata_table_for_all_RDG_samples.tsv")
#geneReadCounts = fread("OUTRIDER_input_table_RDG_counts_for_all_tissues.tsv", select=c('gene_id', sample_subset))
#RDG_subset = RGP_muscle
sample_label = paste("estonia_muscle_", sander_muscle[1], sep="")
RDG_subset = sander_muscle
GTEX_subset = GTEX_muscle_sample_ids

sample_subset = c(RDG_subset, GTEX_subset)

print(sample_subset)
sampleInfo = fread("metadata_table_for_all_RDG_and_GTEX_samples.tsv")
geneReadCounts = fread("OUTRIDER_input_table_RDG_and_GTEX_counts_for_all_tissues.tsv", select=c('gene_id', sample_subset))

sampleInfo = sampleInfo[sampleInfo$sample_id %in% sample_subset]

sampleInfo$read_length = as.character(sampleInfo$read_length)
geneReadCounts = geneReadCounts[!grep('ERCC', geneReadCounts$geneId),]

geneIds = geneReadCounts$gene_id
colsMiusGeneId = colnames(geneReadCounts)[!colnames(geneReadCounts) %in% c('gene_id')]
geneReadCounts = geneReadCounts[,..colsMiusGeneId]
rownames(geneReadCounts) = geneIds

cnts = as.matrix(geneReadCounts)
rownames(cnts) = geneIds
ncol(cnts) 
nrow(cnts)

sampleInfo[,sampleID:=sample_id]
ods <- OutriderDataSet(countData=cnts, colData=sampleInfo)


ods <- estimateSizeFactors(ods)
sortedSizeFactors = sort(sizeFactors(ods))
ggplot(data=NULL, aes(y=sortedSizeFactors, x=1:ncol(ods))) + 
  geom_point(color='blue', size=1) + 
  labs(x='Sample rank', y='Size factors', title="Size factor distribution") + 
  geom_label_repel(aes(label=ifelse(sortedSizeFactors > 1.5, names(sortedSizeFactors), "")), 
                   nudge_x = -35, 
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50') +
  geom_label_repel(aes(label=ifelse(sortedSizeFactors < 0.5, names(sortedSizeFactors), "")), 
                   nudge_x = 35, 
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50') +
  theme_bw()

# from https://i12g-gagneurweb.in.tum.de/public/workshops/RNAseq_ASHG19/input_data/annotations/gencode.v29lift37.annotation.txdb

txdb <- loadDb("/Users/weisburd/project__rnaseq/data/samples/expression/gencode.v26.annotation.txdb")
ods <- filterExpression(ods, gtfFile=txdb, filterGenes=FALSE) #, fpkmCutoff=100)

options(repr.plot.width=4, repr.plot.height=4)
plotFPKM(ods) + theme_bw() + theme(legend.position = 'bottom')

plotExpressedGenes(ods)

ods <- ods[mcols(ods)$passedFilter,]
ods

# Make heatmap figure bigger
options(repr.plot.width=12, repr.plot.height=5)

plotCountCorHeatmap(ods, colGroups=possibleConfounders, normalize=FALSE) # use normalize=FALSE since the data is not yet corrected
#ggsave(file=paste(sample_label, "_corHeatMap", "_uncorrected", ".png", sep=""), width=12, height=8, dpi=150, device="png")


# Heatmap of the gene/sample expression
plotCountGeneSampleHeatmap(ods, colGroups=possibleConfounders, normalized=FALSE)


#ods <- findEncodingDim(ods)
#plotEncDimSearch(ods)

#results_list = list()
original_ods = ods

for(i in 1:4) { 
  start_time = Sys.time()
  print(paste("Start time for q=", 10*i, start_time))
  ods = OUTRIDER(original_ods, verbose=FALSE, iterations=15, q=10*i, # q=35, iterations=10, # only used here 2 iterations to speed up the tutorial
                  BPPARAM=MulticoreParam(4, progressbar=TRUE))  # Takes around 2min
  
  options(repr.plot.width=10, repr.plot.height=10)
  plotCountGeneSampleHeatmap(ods, colGroups=possibleConfounders, normalized=TRUE, main=paste("Count Gene vs Sample Heatmap q=", 10*i, sep=""))
  #ggsave(file=paste(sample_label, "_geneSampleHeatmap_", "q=", 10*i, ".png", sep=""), width=12, height=8, dpi=150, device="png")
  
  plotCountCorHeatmap(ods, colGroups=possibleConfounders, normalize=TRUE, main=paste("Count correlation heatmap q=", 10*i, sep=""))
  #ggsave(file=paste(sample_label, "_corHeatHeatmap_", "q=", 10*i, ".png", sep=""), width=12, height=8, dpi=150, device="png")
  
  saveRDS(ods, paste(sample_label, "_ods__", "q", 10*i, ".RDS", sep=""))
  
  end_time = Sys.time()
  print(paste("End time for q=", 10*i, end_time))
  print(paste("Duration for q=", 10*i, end_time - start_time))
  res <- results(ods)
  head(res)
  print(dim(res))
  
  #results_list = append(results_list, list(ods))
  print("------------------------------------------")
}

options(repr.plot.width=4, repr.plot.height=4)

plotAberrantPerSample(ods)

plotAberrantPerSample(ods, padjCutoff=1, zScoreCutoff=2, yadjust=c(1.3,1.4))

# remove versioning of gene IDs and merge with GRCh37 from annotables
geneIDs <- gsub("\\.[0-9]*(_[0-9]*)?.*$", "", rownames(ods))
map <- merge(data.table(ensgene=geneIDs), grch37, sort=FALSE, all.x=TRUE)[!duplicated(ensgene),]

# set new gene names only if hgnc symbol is present
if(!"ENSG" %in% colnames(mcols(ods))){
  mcols(ods)$ENSG <- geneIDs
  rownames(ods) <- map[,ifelse(is.na(symbol) | symbol == "" | duplicated(symbol), geneIDs, symbol)]
}

res <- results(ods)
head(res)
dim(res)


non_gtex_sample_ids = sampleInfo$sample_id[!grepl('GTEX', sampleInfo$sample_id)]

require(gridExtra)
require(lattice)

plots = lapply(1:length(non_gtex_sample_ids), function(i) {
  print(i)
  sample_id = non_gtex_sample_ids[i]
  print(sample_id)
  volcanoPlotLabels = ifelse(names(ods) %in% res[sampleID == sample_id]$geneID, names(ods), "")
  plotVolcano(ods, sample_id, basePlot=TRUE) +
    geom_label_repel(aes(label=volcanoPlotLabels), force=3, nudge_y = -1, box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50') +
    labs(title=sample_id, x = "", y = "") +
    theme(plot.margin=unit(c(0.5, 0, 0, 0), "cm"))
})

g = grid.arrange(
  grobs=plots, 
  ncol=4, 
  top = textGrob("Estonian Muscle Samples"),
  left = textGrob("-log10(Pvalue)", rot = 90),
  bottom = textGrob("Z-score"),
)

ggsave(file="estonian_muscle_samples_by_pvalue.png", g, width=12, height=8, dpi=150)

options(repr.plot.width=7, repr.plot.height=4)

ggarrange(ncol=2,
          plotExpressionRank(ods, "TIMMDC1", norm=FALSE, basePlot=TRUE),
          plotExpressionRank(ods, "TIMMDC1", norm=TRUE,  basePlot=TRUE))



# make the plot wider
options(repr.plot.width=7, repr.plot.height=4)

# Volcano plot with significance based on p value and on Z score threshold
ggarrange(ncol=2,
          plotVolcano(ods, "MUN_FAM3_PARTIALMDC1A1_05_R1", main="FDR (<0.05) based volcano plot", 
                      base=TRUE, padjCutoff=0.05, zScoreCutoff=0),
          plotVolcano(ods, "MUN_FAM3_PARTIALMDC1A1_05_R1", main="Z score (>3) based volcano plot",
                      base=TRUE, padjCutoff=1, zScoreCutoff=3))


x = list()
for(i in 1:4) { 
  print(i)
  ods = readRDS(paste(sample_label, "_ods__", "q", 10*i, ".RDS", sep=""))
  results(ods)
  g = plotAberrantPerSample(ods, main = paste("Aberrant Genes per Sample q=", 10*i, sep=""))
  x = append(x, list(g))
}

grid.arrange(
  grobs=x, 
  ncol=2, 
  #top = textGrob("Estonian Muscle Samples"),
  #left = textGrob("-log10(Pvalue)", rot = 90),
  #bottom = textGrob("Z-score"),
)

sample_id = "OUN_HK079_001"
sample_id = "OUN_HK018_0047"
sample_id = "OUN_HK137_001"
sample_id = "OUN_HK116_001"
sample_id = "HK072-001_2"

#findEncodingDim(ods, BPPARAM=MulticoreParam(4, progressbar=TRUE))

for(i in 1:4) { 
  ods = readRDS(paste(sample_label, "_ods__", "q", 10*i, ".RDS", sep=""))
  # remove versioning of gene IDs and merge with GRCh37 from annotables
  geneIDs <- gsub("\\.[0-9]*(_[0-9]*)?.*$", "", rownames(ods))
  map <- merge(data.table(ensgene=geneIDs), grch37, sort=FALSE, all.x=TRUE)[!duplicated(ensgene),]
  
  # set new gene names only if hgnc symbol is present
  if(!"ENSG" %in% colnames(mcols(ods))){
    mcols(ods)$ENSG <- geneIDs
    rownames(ods) <- map[,ifelse(is.na(symbol) | symbol == "" | duplicated(symbol), geneIDs, symbol)]
  }
  res <- results(ods, padjCutoff=1)
  print(noquote('-------'))
  print(noquote(paste('q=', 10*i, sep="")))
  print(head(res[res$sampleID == sample_id, c('geneID', 'pValue', 'padjust', 'zScore')][order(pValue),], n=10))
}

res <- results(ods)
head(res)
print(dim(res))
plotAberrantPerSample(ods)

plots = lapply(1:length(non_gtex_sample_ids), function(i) {
  print(i)
  sample_id = non_gtex_sample_ids[i]
  print(sample_id)
  volcanoPlotLabels = ifelse(names(ods) %in% res[sampleID == sample_id]$geneID, names(ods), "")
  plotVolcano(ods, sample_id, basePlot=TRUE) +
    geom_label_repel(aes(label=volcanoPlotLabels), force=3, nudge_y = -1, box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50') +
    labs(title=sample_id, x = "", y = "") +
    theme(plot.margin=unit(c(0.5, 0, 0, 0), "cm"))
})

g = grid.arrange(
  grobs=plots, 
  ncol=4, 
  top = textGrob("Estonian Muscle Samples"),
  left = textGrob("-log10(Pvalue)", rot = 90),
  bottom = textGrob("Z-score"),
)

