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
library(ggplot2)
library(data.table)
library(plotly)
library(stringr)
library(RColorBrewer)
library(ggsci)


setwd("~/project__rnaseq/data/samples/expression/")


# 'rnaseqc_tpm', 'star_counts', 
for (data_type in c('rnaseqc_counts')) {
  if (data_type == 'star_counts') {
    file_paths = list.files(data_type, pattern = "*.tab")
  } else if (data_type == 'rnaseqc_counts') {
    file_paths = list.files(data_type, pattern = "*.gct")
  } else {
    file_paths = list.files(data_type, pattern = "*.gct")
  }
  
  print(paste('# ', data_type))
  df_all = NA
  for (path in file_paths) {
    full_path = paste(data_type, path, sep='/')
    print(paste(path, ": ", full_path))
    sampleName = path
    if (data_type == 'star_counts') {
      sampleName = gsub('.ReadsPerGene.out.tab', '', sampleName)
      df = fread(full_path, col.names = c('geneId', 'counts_unstranded', 'counts_strand1', 'counts_strand2'))
      df = df[geneId != 'N_unmapped' & geneId != 'N_noFeature' & geneId != 'N_multimapping' & geneId != 'N_ambiguous']
      df$counts = df$counts_strand1 + df$counts_strand2
      df = df[ , c('geneId', 'counts'), drop = F]
      colnames(df) = c('geneId', sampleName)
    } else {
      sampleName = gsub('.gene_reads.gct', '', sampleName)
      sampleName = gsub('.gene_tpm.gct', '', sampleName)
      df = fread(full_path)
      df = subset(df, select=-c(Description))
      colnames(df) = c('geneId', sampleName)
    }
    
    if (is.na(df_all)) {
      df_all = df
    } else {
      df_all = merge(df_all, df, by='geneId', all=T)
    }
  }
  
  df = df_all
  
  if (data_type == 'star_counts') {
    df_star = df
  } else if (data_type == 'rnaseqc_counts') {
    df_rnaseqc = df
  } else if (data_type == 'rnaseqc_tpm') {
    df_tpm = df
  }
}

# sample annotations
samples = colnames(df)
samples = samples[samples != 'geneId']
sampleInfo = data.table(samples = samples)
sampleInfo$groups = sampleInfo$samples
sampleInfo$groups = ifelse(grepl('ATYPICAL', sampleInfo$samples), 'PATIENT_ATYPICAL', sampleInfo$groups)
sampleInfo$groups = ifelse(grepl('SIBLING', sampleInfo$samples), 'SIBLING_TOTAL', sampleInfo$groups)
sampleInfo$groups = ifelse(grepl('TOTALMD', sampleInfo$samples), 'TOTAL_CONTROL', sampleInfo$groups)
sampleInfo$groups = ifelse(grepl('PARTIALMD', sampleInfo$samples), 'PARTIAL_CONTROL', sampleInfo$groups)
sampleInfo$groups = ifelse(grepl('FAM1_CTRL', sampleInfo$samples), 'NORMAL_CONTROL', sampleInfo$groups)

sampleInfo$affected = grepl('CONTROL', sampleInfo$groups)
names(sampleInfo) <- sub("samples", "INDIVIDUAL", names(sampleInfo))
sampleInfo

#fwrite(sampleInfo, "sample_annotations.tsv", verbose = T, sep = "\t")
#sampleInfo = fread("sample_annotations.tsv")

df_rnaseqc_counts = as.matrix(df_rnaseqc[, 2:length(colnames(df_rnaseqc))])
rownames(df_rnaseqc_counts) = df_rnaseqc$geneId

ods <- OutriderDataSet(countData=df_rnaseqc_counts, colData=sampleInfo)

ods <- estimateSizeFactors(ods)
ggplot(data=NULL, aes(y=sort(sizeFactors(ods)), x=1:ncol(ods))) + 
  geom_point() + 
  labs(x='Sample rank', y='Size factors', title="Size factor distribution") + 
  theme_bw()

# from https://i12g-gagneurweb.in.tum.de/public/workshops/RNAseq_ASHG19/input_data/annotations/gencode.v29lift37.annotation.txdb


# gencode.v29lift37.annotation.txdb
txdb <- loadDb("/Users/weisburd/project__rnaseq/code/rnaseq_methods/gagneurlab/gencode.v26.annotation.txdb")
txdb <- loadDb("/Users/weisburd/project__rnaseq/code/rnaseq_methods/gagneurlab/gencode.v29lift37.annotation.txdb")
txdb <- loadDb("/Users/weisburd/project__rnaseq/code/rnaseq_methods/gagneurlab/gencode.v32.annotation.txdb")
ods <- filterExpression(ods, gtfFile=txdb, filterGenes=FALSE, fpkmCutoff=100)

plotFPKM(ods) + theme_bw() + theme(legend.position = 'bottom')

ods <- ods[mcols(ods)[,"passedFilter"]]
ods

# Make heatmap figure bigger
options(repr.plot.width=6, repr.plot.height=5)

# use normalize=FALSE since the data is not yet corrected
plotCountCorHeatmap(ods, colGroups=c("groups"), normalize=FALSE)


ods <- OUTRIDER(ods, q=15, verbose=TRUE, iterations=2, # only used here 2 iterations to speed up the tutorial
                BPPARAM=MulticoreParam(4, progressbar=TRUE))  # Takes around 2min




plotCountCorHeatmap(ods, rowGroups="groups",normalize=TRUE)


options(repr.plot.width=4, repr.plot.height=4)

plotAberrantPerSample(ods)

plotAberrantPerSample(ods, padjCutoff=1, zScoreCutoff=2, yadjust=c(1.3,1.4))

res <- results(ods)
head(res)
dim(res)


plotVolcano(ods, "MUN_FAM3_PARTIALMDC1A1_05_R1", base=TRUE)


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


# PCA
df_tpm_with_geneId = df_tpm


df_tpm = subset(df_tpm, select=-c(geneId))
df_tpm = df_tpm[rowSums(df_tpm) > 0,]
PCA = prcomp(log2(t(as.matrix(df_tpm)) + 1), center=T, retx = T)
PCA_df = data.frame(PCA$x)
PCA_df$samples = rownames(PCA_df)

PCA_df = merge(PCA_df, sampleInfo, by='samples', all=T)

ggplot(PCA_df, aes(x=PC1, y=PC2, col=groups) ) + geom_point(size=1.5)
ggplot(PCA_df, aes(x=PC1, y=PC3, col=groups) ) + geom_point(size=1.5)
ggplot(PCA_df, aes(x=PC2, y=PC3, col=groups) ) + geom_point(size=1.5)

plot_ly(PCA_df, x=~PC1, y=~PC2, z=~PC3, text=~samples, color=~groups)
plot_ly(PCA_df, x=~PC1, y=~PC2, z=~PC4, color=~groups)

M = log2(t(as.matrix(df_tpm)) + 1)
d = dist(scale(M), method = "euclidean")
hc = hclust(d, method = "ward.D")  # complete, ward.D, centroid, etc.
#hc = hclust(d, method = "complete")  # complete, ward.D, centroid, etc.

plot(hc, cex = 0.6, hang = -1)
rect.hclust(hc, k = 8, border = 1:7)

# GTEx
gtex_sampleInfo = fread('~/project__rnaseq/data/GTEx_v8_metadata/GTEx_Analysis_v8_Annotations_SampleAttributesDS.tsv')

gtex_individualInfo = fread('~/project__rnaseq/data/GTEx_v8_metadata/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.tsv')
gtex_sampleInfo$SUBJID = unlist(lapply(gtex_sampleInfo$SAMPID, function(x) { paste(str_split_fixed(x, '-', 3)[1:2], collapse='-') }))

gtex_sampleInfo = merge(gtex_sampleInfo, gtex_individualInfo, by='SUBJID', all.x=T)
gtex_sampleInfo$SEX = ifelse(gtex_sampleInfo$SEX == 1, 'MALE', ifelse(gtex_sampleInfo$SEX == 2, 'FEMALE', NA))

broad_tissues = c('Blood', 'Adipose Tissue', 'Muscle', 'Blood Vessel', 'Heart', 'Breast', 'Skin', 'Lung', 'Esophagus', 'Nerve')
gtex_samples = gtex_sampleInfo[(gtex_sampleInfo$SMTS %in% broad_tissues) & gtex_sampleInfo$SMRIN >= 5.7]$SAMPID

narrow_tissues = c('Blood', 'Adipose Tissue', 'Muscle', 'Skin',  'Blood Vessel', 'Nerve')
gtex_samples = gtex_sampleInfo[(gtex_sampleInfo$SMTS %in% narrow_tissues) & gtex_sampleInfo$SMRIN >= 5.7]$SAMPID

#focused_tissue_subset = c('Adipose - Subcutaneous', 'Muscle - Skeletal')
#gtex_samples = gtex_sampleInfo[(gtex_sampleInfo$SMTSD %in% focused_tissue_subset) & gtex_sampleInfo$SMRIN > 5.7]$SAMPID

print(paste(length(gtex_samples), 'GTEx samples kept'))

#gtex_samples = gtex_samples[!gtex_samples %in% c('GTEX-5E5C-1626-SM-2XCE3')]  # exclude sample mismatches

gtex_tpm = fread(
  '~/project__rnaseq/data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz',
  select=gtex_samples)
#drop=c('Names', 'Description'))

gtex_tpm_names = fread(
  '~/project__rnaseq/data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz',
  select=c('Name'))


# GTEx PCA
gtex_tpm_with_geneId = gtex_tpm
gtex_tpm_with_geneId$geneId = gtex_tpm_names$Name


gtex_tpm = gtex_tpm[rowSums(gtex_tpm) > 0,]
gtex_PCA = prcomp(log2(t(as.matrix(gtex_tpm)) + 1), center=T, retx = T)
gtex_PCA_df = data.frame(gtex_PCA$x)
gtex_PCA_df$SAMPID = row.names(gtex_PCA_df)

gtex_PCA_df = merge(gtex_PCA_df[,c('SAMPID', 'PC1', 'PC2', 'PC3')], gtex_sampleInfo[,c('SAMPID','SMTS', 'SMTSD', 'SMAFRZE', 'SMRIN')], by='SAMPID', all.x=T)

ggplot(gtex_PCA_df, aes(x=PC1, y=PC2, col=SMTSD) ) + geom_point(size=1.5)

plot_ly(gtex_PCA_df, x=~PC1, y=~PC2, z=~PC3, color=~SMTSD, text=~SAMPID)


# merge GTEx + samples
merged_tpm = merge(df_tpm_with_geneId, gtex_tpm_with_geneId, by="geneId")
colsMiusGeneId = colnames(merged_tpm)[!colnames(merged_tpm) %in% c('geneId')]
merged_tpm = merged_tpm[,..colsMiusGeneId]


merged_tpm = merged_tpm[rowSums(merged_tpm) > 0,]
merged_PCA = prcomp(log2(t(as.matrix(merged_tpm)) + 1), center=T, retx = T)
merged_PCA_df = data.frame(merged_PCA$x)
merged_PCA_df$SAMPID = row.names(merged_PCA_df)

merged_PCA_df = merge(merged_PCA_df[,c('SAMPID', 'PC1', 'PC2', 'PC3')], sampleInfo, by.x='SAMPID', by.y='samples', all.x=T)
merged_PCA_df = merge(merged_PCA_df[,c('SAMPID', 'PC1', 'PC2', 'PC3', 'groups')], gtex_sampleInfo[,c('SAMPID','SMTS', 'SMTSD', 'SMAFRZE', 'SMRIN', 'SEX', 'AGE', 'DTHHRDY')], by='SAMPID', all.x=T)
merged_PCA_df$label = ifelse(is.na(merged_PCA_df$groups), merged_PCA_df$SMTSD, merged_PCA_df$groups)

ggplot(merged_PCA_df, aes(x=PC1, y=PC2, col=SEX) ) + geom_point(size=1.5)  + scale_color_lancet()

plot_ly(
  merged_PCA_df,
  x=~PC1, y=~PC2, z=~PC3, color=~DTHHRDY, text=~SAMPID,
  marker = list(size = 5))

plot_ly(
  merged_PCA_df,
  x=~PC1, y=~PC2, z=~PC3, color=~label, text=~SAMPID,
  colors='Accent', marker = list(size = 5))


nearby_GTEx_samples = c(
  'GTEX-13YAN-0526-SM-5O9BE',
  'GTEX-139T4-0326-SM-5K7XN',
  'GTEX-183WM-0426-SM-731B1',
  'GTEX-1LVAO-0526-SM-EWRNS',
  'GTEX-WOFL-0626-SM-3MJG3',
  'GTEX-13X6I-0526-SM-5QGP8',
  'GTEX-1R9K5-0326-SM-EWRNG',
  'GTEX-ZZPT-0626-SM-5GZXT',
  'GTEX-144FL-0326-SM-5O99O',
  'GTEX-14ASI-0526-SM-5QGQP',
  
  'GTEX-1JMPZ-0126-SM-CNPOZ',
  'GTEX-X88G-0326-SM-47JZ4',
  'GTEX-VUSG-2626-SM-4KKZI',
  'GTEX-WWTW-0626-SM-4MVNS',
  'GTEX-1CB4H-0426-SM-7MXTL',
  'GTEX-1ICLZ-0626-SM-ARL83',
  'GTEX-OHPN-2726-SM-2I5H4',
  'GTEX-1F5PK-2626-SM-7RHH6',
  'GTEX-132AR-1026-SM-5PNVL',
  
  'GTEX-QDVN-2426-SM-2S1Q4',
  'GTEX-17EVQ-1226-SM-7LG59',
  'GTEX-XMD1-0626-SM-EZ6MC',
  'GTEX-1F88E-1826-SM-9MQLV',
  'GTEX-1J8EW-0326-SM-CYPSX',
  'GTEX-1KANB-0126-SM-E6CJT',
  'GTEX-13SLW-0326-SM-5RQK5',
  'GTEX-SSA3-0326-SM-32QPS',
  'GTEX-13IVO-0626-SM-5LZYJ',
  'GTEX-13113-5019-SM-7EPH2',
  'GTEX-1IGQW-0526-SM-CKZO3',
  'GTEX-YBZK-0326-SM-59HLN',
  'GTEX-1H2FU-0526-SM-ACKXO',
  'GTEX-1J1OQ-0426-SM-A9SLI'
)



# subset

# merge GTEx + samples
merged_nearby_tpm = merge(df_tpm_with_geneId, subset(gtex_tpm_with_geneId, select=c('geneId', nearby_GTEx_samples)), by="geneId")
colsMiusGeneId = colnames(merged_nearby_tpm)[!colnames(merged_nearby_tpm) %in% c('geneId')]
merged_nearby_tpm = merged_nearby_tpm[,..colsMiusGeneId]


merged_nearby_tpm = merged_nearby_tpm[rowSums(merged_nearby_tpm) > 0,]
merged_nearby_PCA = prcomp(log2(t(as.matrix(merged_nearby_tpm)) + 1), center=T, retx = T)
merged_nearby_PCA_df = data.frame(merged_nearby_PCA$x)
merged_nearby_PCA_df$SAMPID = row.names(merged_nearby_PCA_df)

merged_nearby_PCA_df = merge(merged_nearby_PCA_df[,c('SAMPID', 'PC1', 'PC2', 'PC3')], sampleInfo, by.x='SAMPID', by.y='samples', all.x=T)
merged_nearby_PCA_df = merge(merged_nearby_PCA_df[,c('SAMPID', 'PC1', 'PC2', 'PC3', 'groups')], gtex_sampleInfo[,c('SAMPID','SMTS', 'SMTSD', 'SMAFRZE', 'SMRIN', 'SEX', 'AGE', 'DTHHRDY')], by='SAMPID', all.x=T)
merged_nearby_PCA_df$label = ifelse(is.na(merged_nearby_PCA_df$groups), merged_nearby_PCA_df$SMTSD, merged_nearby_PCA_df$groups)

ggplot(merged_nearby_PCA_df, aes(x=PC1, y=PC2, col=SEX) ) + geom_point(size=1.5)  + scale_color_lancet()


merged_nearby_PCA_df$SEX2 = ifelse(is.na(merged_nearby_PCA_df$SEX), 'UNKNOWN', merged_nearby_PCA_df$SEX)

plot_ly(
  merged_nearby_PCA_df,
  x=~PC1, y=~PC2, z=~PC3, color=~DTHHRDY, text=~SAMPID,
  marker = list(size = 5))
plot_ly(
  merged_nearby_PCA_df,
  x=~PC1, y=~PC2, z=~PC3, color=~label, text=~SAMPID,
  colors='Set2', marker = list(size = 5))

plot_ly(
  merged_nearby_PCA_df,
  x=~PC1, y=~PC2, z=~PC3, color=~SEX2, text=~SAMPID,
  colors='Accent', marker = list(size = 5))

M = log2(t(as.matrix(merged_nearby_tpm)) + 1)
d = dist(scale(M), method = "euclidean")
hc = hclust(d, method = "ward.D")  # complete, ward.D, centroid, etc.
#hc = hclust(d, method = "complete")  # complete, ward.D, centroid, etc.

plot(hc, cex = 0.6, hang = -1)
rect.hclust(hc, k = 5, border = 1:4)

