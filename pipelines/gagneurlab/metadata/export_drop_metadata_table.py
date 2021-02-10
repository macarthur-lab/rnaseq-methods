import argparse
import datetime
import glob
import gzip
import hashlib
import logging
import os
import pandas as pd
import sys

from sample_metadata.rnaseq_metadata_utils import get_joined_metadata_df, get_gtex_rnaseq_sample_metadata_df

logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s', level=logging.INFO)
logger = logging.getLogger(__name__)

RDG_GENE_READS_PATH = os.path.expanduser("~/project__rnaseq/data/samples/expression/rnaseqc_counts/*gene_reads.gct.gz")
GTEX_GENE_READS_PATH = os.path.expanduser('~/project__rnaseq/data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz')


def get_sample_set_label(sample_ids):
    byte_string = ", ".join(sorted(sample_ids)).encode()
    h = hashlib.md5(byte_string).hexdigest().upper()
    return f"{len(sample_ids)}_samples_{h[:10]}"


#  RNA_ID
#  RNA_BAM_FILE
#  DNA_VCF_FILE
#  DNA_ID
#  DROP_GROUP - example: "import_exp", "outrider,mae", "fraser"
#  PAIRED_END - example: "True"
#  COUNT_MODE - example: "IntersectionStrict"
#  COUNT_OVERLAPS - example: "True"
#  STRAND - example: "no"
#  HPO_TERMS
#  GENE_COUNTS_FILE
#  ANNOTATION - example: "v29"

def transfer_metadata_columns_from_df(samples_df, source_df):
    df = pd.DataFrame()
    df.loc[:, 'RNA_ID'] = source_df['sample_id']
    df.loc[:, 'DNA_ID'] = source_df['sample_id']
    #df.loc[source_df.sample_id, 'batch'] = source_df['star_pipeline_batch']
    #df.loc[source_df.sample_id, 'batch_detail'] = source_df['batch_date_from_hg19_bam_header']
    df.loc[source_df.sample_id, 'RNA_BAM_FILE'] = source_df['star_bam']
    df.loc[source_df.sample_id, 'RNA_BAI_FILE'] = source_df['star_bai']
    df.loc[source_df.sample_id, 'DNA_VCF_FILE'] = source_df['grch38_vcf']
    df.loc[source_df.sample_id, 'PAIRED_END'] = True
    df.loc[source_df.sample_id, 'STRAND'] = source_df['stranded? (rnaseqc)']  #  "yes", "no"
    df.loc[source_df.sample_id, 'HPO_TERMS'] = ""
    df.loc[source_df.sample_id, 'COUNT_MODE'] = "IntersectionStrict"
    df.loc[source_df.sample_id, 'COUNT_OVERLAPS'] = True
    df.loc[source_df.sample_id, 'ANNOTATION'] = "v34"
    df.loc[source_df.sample_id, 'GENE_COUNTS_FILE'] = ""

    #df.loc[source_df.sample_id, 'star_SJ_out_tab'] = source_df['star_SJ_out_tab']
    #df.loc[source_df.sample_id, 'star_reads_per_gene_tab'] = source_df['star_reads_per_gene_tab']
    #df.loc[source_df.sample_id, 'rnaseqc_gene_reads'] = source_df['rnaseqc_gene_reads']
    #df.loc[source_df.sample_id, 'rnaseqc_gene_tpm'] = source_df['rnaseqc_gene_tpm']

    #df.loc[source_df.sample_id, 'output_dir'] = source_df['star_pipeline_batch'].apply(lambda batch_name: f"gs://macarthurlab-rnaseq/{batch_name}/fraser_count_rna/")
    #df.loc[source_df.sample_id, 'project'] = source_df['proj (seqr)']
    #df.loc[source_df.sample_id, 'sex'] = source_df['imputed sex']
    #df.loc[source_df.sample_id, 'age'] = source_df['Age at muscle biopsy (Beryl:Supp.)'].apply(lambda x: ('' if not x or not x.strip() else ('0' if (x and (x == 'At birth' or 'm' in x)) else x.replace('y', ''))))
    #df.loc[source_df.sample_id, 'cause_of_death'] = None
    #df.loc[source_df.sample_id, 'ancestry'] = None
    #df.loc[source_df.sample_id, 'tissue'] = None
    #df.loc[source_df.sample_id, 'tissue_detail'] = None
    #df.loc[source_df.sample_id, 'read_length'] = source_df['read length (rnaseqc)']
    #df.loc[source_df.sample_id, 'stranded'] = source_df['stranded? (rnaseqc)']
    #df.loc[source_df.sample_id, 'RIN'] = None

    return pd.concat([samples_df, df], axis="rows", sort=True)


def transfer_metadata_columns_from_GTEx_df(samples_df, source_df, batch_name):
    print("Adding")
    print(source_df[['SAMPID', 'SMRIN']])

    df = pd.DataFrame()
    df.loc[:, 'RNA_ID'] = source_df['SAMPID']
    df.loc[:, 'DNA_ID'] = source_df['SUBJID']
    #df.loc[source_df.SAMPID, 'batch'] = 'GTEx_v8'
    df.loc[source_df.SAMPID, 'RNA_BAM_FILE'] = source_df['rnaseq_bam']
    df.loc[source_df.SAMPID, 'RNA_BAI_FILE'] = source_df['rnaseq_bai']
    #df.loc[source_df.SAMPID, 'DNA_VCF_FILE'] = source_df['wes_vcf']
    df.loc[source_df.SAMPID, 'DNA_VCF_FILE'] = source_df['wgs_vcf']
    df.loc[source_df.SAMPID, 'PAIRED_END'] = True
    df.loc[source_df.SAMPID, 'STRAND'] = "no"
    df.loc[source_df.SAMPID, 'HPO_TERMS'] = ""
    df.loc[source_df.SAMPID, 'COUNT_MODE'] = "IntersectionStrict"
    df.loc[source_df.SAMPID, 'COUNT_OVERLAPS'] = True
    df.loc[source_df.SAMPID, 'ANNOTATION'] = "v34"
    df.loc[source_df.SAMPID, 'GENE_COUNTS_FILE'] = ""

    #df.loc[source_df.SAMPID, 'batch_detail'] = "GTEx_v8 " + source_df['SMNABTCH']
    #df.loc[source_df.SAMPID, 'output_dir'] = f"gs://macarthurlab-rnaseq/gtex_v8/fraser_count_rna/"
    #df.loc[source_df.SAMPID, 'project'] = "gtex_v8"
    #df.loc[source_df.SAMPID, 'sex'] = source_df['SEX']
    #df.loc[source_df.SAMPID, 'age'] = source_df['AGE']
    #df.loc[source_df.SAMPID, 'cause_of_death'] = source_df['DTHHRDY']
    #df.loc[source_df.SAMPID, 'ancestry'] = None
    #df.loc[source_df.SAMPID, 'tissue'] = source_df['SMTS']
    #df.loc[source_df.SAMPID, 'tissue_detail'] = source_df['SMTSD']
    #df.loc[source_df.SAMPID, 'read_length'] = source_df['SMRDLGTH']
    #df.loc[source_df.SAMPID, 'RIN'] = source_df['SMRIN']

    return pd.concat([samples_df, df], axis="rows")


TISSUE_NAME_TO_SMTD = {
    "muscle": "Muscle - Skeletal",
    "fibroblasts": "Cells - Cultured fibroblasts",
    "whole_blood": "Whole Blood",
    "lymphocytes": "Cells - EBV-transformed lymphocytes",
}


def main():
    # get metadata table and filter out excluded samples
    rnaseq_sample_metadata_df = get_joined_metadata_df()
    rnaseq_sample_metadata_df = rnaseq_sample_metadata_df[~rnaseq_sample_metadata_df["analysis batch"].str.strip().isin(["", "x"])]

    gtex_rnaseq_sample_metadata_df = get_gtex_rnaseq_sample_metadata_df()

    all_rdg_samples_df = None

    all_rdg_and_gtex_samples_df = None
    all_rdg_and_gtex_counts_df = None

    for tissue_name in ["whole_blood", "lymphocytes", "fibroblasts", "muscle"]:
        logging.info("----------------")
        logging.info(f"{tissue_name}:")
        samples_df = pd.DataFrame()
        SMTD_value = TISSUE_NAME_TO_SMTD[tissue_name]

        # process RDG samples
        tgg_df = rnaseq_sample_metadata_df[rnaseq_sample_metadata_df['imputed tissue'] == tissue_name]
        logging.info(f"Got {len(tgg_df)} TGG samples for {tissue_name}")
        samples_df = transfer_metadata_columns_from_df(samples_df, tgg_df)
        drop_groups = f"{tissue_name},{tissue_name}_with_GTEx" if len(tgg_df) > 30 else f"{tissue_name}_with_GTEx"
        samples_df.loc[tgg_df.sample_id, 'DROP_GROUP'] = drop_groups

        #samples_df['tissue_detail'] = SMTD_value
        if all_rdg_samples_df is None:
            all_rdg_samples_df = samples_df
        else:
            all_rdg_samples_df = pd.concat([all_rdg_samples_df, samples_df])

        # process GTEx samples
        gtex_df = gtex_rnaseq_sample_metadata_df[gtex_rnaseq_sample_metadata_df.SMTSD == SMTD_value]
        gtex_df = gtex_df[gtex_df["wgs_vcf"].str.len() > 1]
        gtex_df = gtex_df.sort_values(by='SMRIN', ascending=False)
        logging.info(f"Got {len(gtex_df)} GTEx samples for {tissue_name}")

        samples_df = transfer_metadata_columns_from_GTEx_df(samples_df, gtex_df[:100], tissue_name)
        samples_df.loc[gtex_df[:100].SAMPID, 'DROP_GROUP'] = f"{tissue_name}_with_GTEx"

        if all_rdg_and_gtex_samples_df is None:
            all_rdg_and_gtex_samples_df = samples_df
        else:
            all_rdg_and_gtex_samples_df = pd.concat([all_rdg_and_gtex_samples_df, samples_df])

        #tsv_output_path = f"metadata_table_for_{tissue_name}.tsv"
        #samples_df.to_csv(tsv_output_path, sep="\t", index=False)
        #print(f"Wrote {len(samples_df)} samples to {tsv_output_path}")

    # create separate groups for MAE analysis, including only samples that have a VCF
    all_rdg_and_gtex_samples_df.loc[
        all_rdg_and_gtex_samples_df["DNA_VCF_FILE"].str.len() > 1,
        'DROP_GROUP'] = all_rdg_and_gtex_samples_df['DROP_GROUP'].apply(
            lambda groups: groups+","+",".join([f"{g}_MAE" for g in groups.split(",")]))

    #  RNA_ID
    #  RNA_BAM_FILE
    #  DNA_VCF_FILE
    #  DNA_ID
    #  DROP_GROUP - example: "import_exp", "outrider,mae", "fraser"
    #  PAIRED_END - example: "True"
    #  COUNT_MODE - example: "IntersectionStrict"
    #  COUNT_OVERLAPS - example: "True"
    #  STRAND - example: "no"
    #  HPO_TERMS
    #  GENE_COUNTS_FILE
    #  ANNOTATION - example: "v29"
    all_rdg_and_gtex_samples_df = all_rdg_and_gtex_samples_df[[
        "RNA_ID", "DNA_ID", "RNA_BAM_FILE", "DNA_VCF_FILE", "DROP_GROUP", "PAIRED_END",
        "COUNT_MODE", "COUNT_OVERLAPS", "STRAND", "HPO_TERMS", "GENE_COUNTS_FILE", "ANNOTATION",
    ]]
    print("Output columns:")
    print(all_rdg_and_gtex_samples_df.columns)

    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    tsv_output_path = f"drop_metadata_table__{timestamp}.tsv"
    all_rdg_and_gtex_samples_df.to_csv(tsv_output_path, sep="\t", index=False)
    print(f"Wrote {len(all_rdg_and_gtex_samples_df)} samples to {tsv_output_path}")


    #os.system(f"gsutil -m cp {tsv_output_path} gs://seqr-bw/project__rnaseq/")


if __name__ == "__main__":
    main()
