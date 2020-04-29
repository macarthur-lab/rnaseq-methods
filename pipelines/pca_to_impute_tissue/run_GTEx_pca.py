import argparse
import gzip
import logging
import os
import pandas as pd
import pprint

logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s')
logger = logging.getLogger()
logger.setLevel(logging.INFO)

GTEX_SAMPLE_INFO = os.path.expanduser('~/project__rnaseq/data/GTEx_v8_metadata/GTEx_Analysis_v8_Annotations_SampleAttributesDS.tsv')
GTEX_INDIVIDUAL_INFO = os.path.expanduser('~/project__rnaseq/data/GTEx_v8_metadata/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.tsv')
GTEX_TPM_PATH = os.path.expanduser("~/project__rnaseq/data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz")

def parse_args():
    p = argparse.ArgumentParser()
    #p.add_argument("--source-namespace", default="macarthurlab-rnaseq-terra")
    #p.add_argument("batch_name")
    args = p.parse_args()

    return args


def run(cmd):
    logger.info(cmd)
    os.system(cmd)


def gsutil_cp(source, dest):
    run("gsutil -m cp -n %s  %s" % (source, dest))


def main():
    args = parse_args()
    logger.info("Args:\n" + pprint.pformat(args.__dict__))

    df_sampleInfo = pd.read_csv(GTEX_SAMPLE_INFO, sep="\t")
    df_individualInfo = pd.read_csv(GTEX_INDIVIDUAL_INFO, sep="\t")
    with gzip.open(GTEX_TPM_PATH, "rt") as f:
        next(f)  # skip header
        next(f)
        df = pd.read_csv(f, sep="\t")

    df_sampleInfo["SUBJID"] = df_sampleInfo["SAMPID"].apply(lambda x: "-".join(x.split('-')[1:3]))
    print(df_sampleInfo["SUBJID"])

    logger.info("Done")


if __name__ == "__main__":
    main()


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
