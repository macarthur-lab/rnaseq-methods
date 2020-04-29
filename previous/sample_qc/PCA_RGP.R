library(dplyr)
library(purrr)
library(ggplot2)
library(data.table)
library(plotly)
library(stringr)
library(RColorBrewer)
library(ggsci)


setwd("~/project__rnaseq/data/samples/expression/")
data_type = "rnaseqc_tpm"
file_paths = list.files(data_type, pattern = "RGP_.*.gct")
  
print(paste('# ', data_type))
df_all = NA
for (path in file_paths) {
  full_path = paste(data_type, path, sep='/')
  print(paste(path, ": ", full_path))
  sampleName = path
  sampleName = gsub('.gene_tpm.gct', '', sampleName)
  df = fread(full_path)
  df = subset(df, select=-c(Description))
  colnames(df) = c('geneId', sampleName)

  if (is.na(df_all)) {
    df_all = df
  } else {
    df_all = merge(df_all, df, by='geneId', all=T)
  }
}

df = df_all
df_tpm = df



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
#gtex_tpm_with_geneId = gtex_tpm
#gtex_tpm_with_geneId$geneId = gtex_tpm_names$Name


#gtex_tpm = gtex_tpm[rowSums(gtex_tpm) > 0,]
#gtex_PCA = prcomp(log2(t(as.matrix(gtex_tpm)) + 1), center=T, retx = T)
#gtex_PCA_df = data.frame(gtex_PCA$x)
#gtex_PCA_df$SAMPID = row.names(gtex_PCA_df)

#gtex_PCA_df = merge(gtex_PCA_df[,c('SAMPID', 'PC1', 'PC2', 'PC3')], gtex_sampleInfo[,c('SAMPID','SMTS', 'SMTSD', 'SMAFRZE', 'SMRIN')], by='SAMPID', all.x=T)

#ggplot(gtex_PCA_df, aes(x=PC1, y=PC2, col=SMTSD) ) + geom_point(size=1.5)

#plot_ly(gtex_PCA_df, x=~PC1, y=~PC2, z=~PC3, color=~SMTSD, text=~SAMPID)


# merge GTEx + samples
df_tpm_with_geneId = df_tpm
gtex_tpm_with_geneId = gtex_tpm
gtex_tpm_with_geneId$geneId = gtex_tpm_names$Name

merged_tpm = merge(df_tpm_with_geneId, gtex_tpm_with_geneId, by="geneId")
colsMiusGeneId = colnames(merged_tpm)[!colnames(merged_tpm) %in% c('geneId')]
merged_tpm = merged_tpm[,..colsMiusGeneId]


merged_tpm = merged_tpm[rowSums(merged_tpm) > 0,]
merged_PCA = prcomp(log2(t(as.matrix(merged_tpm)) + 1), center=T, retx = T)
merged_PCA_df = data.frame(merged_PCA$x)
merged_PCA_df$SAMPID = row.names(merged_PCA_df)

merged_PCA_df = merge(merged_PCA_df[,c('SAMPID', 'PC1', 'PC2', 'PC3')], gtex_sampleInfo[,c('SAMPID','SMTS', 'SMTSD', 'SMAFRZE', 'SMRIN', 'SEX', 'AGE', 'DTHHRDY')], by='SAMPID', all.x=T)

merged_PCA_df$SMTSD = ifelse(grepl('RGP_248_3', merged_PCA_df$SAMPID), ' RGP_248_3', merged_PCA_df$SMTSD)
merged_PCA_df$SMTSD = ifelse(grepl('RGP_273_3', merged_PCA_df$SAMPID), ' RGP_273_3', merged_PCA_df$SMTSD)
merged_PCA_df$SMTS = ifelse(grepl('RGP_248_3', merged_PCA_df$SAMPID), ' RGP_248_3', merged_PCA_df$SMTS)
merged_PCA_df$SMTS = ifelse(grepl('RGP_273_3', merged_PCA_df$SAMPID), ' RGP_273_3', merged_PCA_df$SMTS)
merged_PCA_df$label = ifelse(grepl('RGP_248_3', merged_PCA_df$SAMPID), ' RGP_248_3', merged_PCA_df$SMTSD)
merged_PCA_df$label = ifelse(grepl('RGP_273_3', merged_PCA_df$SAMPID), ' RGP_273_3', merged_PCA_df$SMTSD)

merged_PCA_df_subset = merged_PCA_df[merged_PCA_df$SMTS %in% c('Blood Vessel',  'Nerve',  'Blood', 'Muscle', 'Skin',  'Adipose Tissue', ' RGP_248_3', ' RGP_273_3'),]

ggplot(merged_PCA_df_subset, aes(x=PC1, y=PC2, col=SMTSD) ) + geom_point(size=1.5) # + scale_color_lancet()

plot_ly(
  merged_PCA_df,
  x=~PC1, y=~PC2, z=~PC3, color=~SMTS, text=~SAMPID,
  colors = c('purple', 'purple', 'gold', '#AA3333', '#00BB00', 'pink', '#AAAAAA', '#00BBCC'),
  marker = list(size = 6))

]plot_ly(
  merged_PCA_df,
  x=~PC1, y=~PC2, z=~PC3, color=~label, text=~SAMPID,
  colors='Accent', marker = list(size = 3))


narrow_tissues = c('Blood', 'Adipose Tissue', 'Muscle', 'Skin',  'Blood Vessel', 'Nerve')
gtex_samples = gtex_sampleInfo[(gtex_sampleInfo$SMTSD %in% narrow_tissues) & gtex_sampleInfo$SMRIN >= 5.7]$SAMPID

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



##--------------------------
# gene list

gtex_sampleInfo = fread('~/project__rnaseq/data/GTEx_v8_metadata/GTEx_Analysis_v8_Annotations_SampleAttributesDS.tsv')
gtex_individualInfo = fread('~/project__rnaseq/data/GTEx_v8_metadata/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.tsv')
gtex_sampleInfo$SUBJID = unlist(lapply(gtex_sampleInfo$SAMPID, function(x) { paste(str_split_fixed(x, '-', 3)[1:2], collapse='-') }))
gtex_sampleInfo = merge(gtex_sampleInfo, gtex_individualInfo, by='SUBJID', all.x=T)
gtex_sampleInfo$SEX = ifelse(gtex_sampleInfo$SEX == 1, 'MALE', ifelse(gtex_sampleInfo$SEX == 2, 'FEMALE', NA))

narrow_tissues = c('Muscle')
gtex_samples = gtex_sampleInfo[(gtex_sampleInfo$SMTS %in% narrow_tissues) & gtex_sampleInfo$SMRIN >= 5.7]$SAMPID

focused_tissue_subset = c('Muscle - Skeletal')
gtex_samples = gtex_sampleInfo[(gtex_sampleInfo$SMTSD %in% focused_tissue_subset) & gtex_sampleInfo$SMRIN > 5.7]$SAMPID

print(paste(length(gtex_samples), 'GTEx samples kept'))

#gtex_samples = gtex_samples[!gtex_samples %in% c('GTEX-5E5C-1626-SM-2XCE3')]  # exclude sample mismatches

gtex_tpm = fread(
  '~/project__rnaseq/data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz',
  select=gtex_samples)
#drop=c('Names', 'Description'))

gtex_tpm_names = fread(
  '~/project__rnaseq/data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz',
  select=c('Name'))


merged_PCA_df_subset2 = merged_PCA_df[merged_PCA_df$SMTS %in% c('Muscle'),]

