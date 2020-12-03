import argparse
import datetime
import glob
import gzip
import hashlib
import logging
import os
import pandas as pd
import sys

from sample_metadata.rnaseq_metadata_utils import get_joined_metadata_df, get_gtex_rnaseq_sample_metadata_df, ANALYSIS_BATCHES

logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s', level=logging.INFO)
logger = logging.getLogger(__name__)

RDG_GENE_READS_PATH = os.path.expanduser("~/project__rnaseq/data/samples/expression/rnaseqc_counts/*gene_reads.gct.gz")
GTEX_GENE_READS_PATH = os.path.expanduser('~/project__rnaseq/data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz')

EXPORT_100_MALE_AND_100_FEMALE_GTEX_SAMPLES = False

def get_sample_set_label(sample_ids):
    byte_string = ", ".join(sorted(sample_ids)).encode()
    h = hashlib.md5(byte_string).hexdigest().upper()
    return f"{len(sample_ids)}_samples_{h[:10]}"


def transfer_metadata_columns_from_df(samples_df, source_df):
    df = pd.DataFrame()
    df.loc[:, 'sample_id'] = source_df['sample_id']
    df.loc[source_df.sample_id, 'batch'] = source_df['star_pipeline_batch']
    df.loc[source_df.sample_id, 'batch_detail'] = source_df['batch_date_from_hg19_bam_header']
    df.loc[source_df.sample_id, 'bam_path'] = source_df['star_bam']
    df.loc[source_df.sample_id, 'bai_path'] = source_df['star_bai']
    df.loc[source_df.sample_id, 'star_SJ_out_tab'] = source_df['star_SJ_out_tab']
    df.loc[source_df.sample_id, 'star_reads_per_gene_tab'] = source_df['star_reads_per_gene_tab']
    df.loc[source_df.sample_id, 'rnaseqc_gene_reads'] = source_df['rnaseqc_gene_reads']
    df.loc[source_df.sample_id, 'rnaseqc_gene_tpm'] = source_df['rnaseqc_gene_tpm']

    #df.loc[source_df.sample_id, 'output_dir'] = source_df['star_pipeline_batch'].apply(lambda batch_name: f"gs://macarthurlab-rnaseq/{batch_name}/fraser_count_rna/")
    df.loc[source_df.sample_id, 'project'] = source_df['proj (seqr)']
    df.loc[source_df.sample_id, 'sex'] = source_df['imputed sex']
    df.loc[source_df.sample_id, 'age'] = source_df['Age at muscle biopsy (Beryl:Supp.)'].apply(lambda x: ('' if not x or not x.strip() else ('0' if (x and (x == 'At birth' or 'm' in x)) else x.replace('y', ''))))
    df.loc[source_df.sample_id, 'cause_of_death'] = None
    df.loc[source_df.sample_id, 'ancestry'] = None
    df.loc[source_df.sample_id, 'tissue'] = None
    df.loc[source_df.sample_id, 'tissue_detail'] = None
    df.loc[source_df.sample_id, 'read_length'] = source_df['read length (rnaseqc)']
    df.loc[source_df.sample_id, 'stranded'] = source_df['stranded? (rnaseqc)']
    df.loc[source_df.sample_id, 'RIN'] = None

    return pd.concat([samples_df, df], axis="rows", sort=True)


def transfer_metadata_columns_from_GTEx_df(samples_df, source_df, batch_name):
    print("Adding")
    print(source_df[['SAMPID', 'SMRIN']])

    df = pd.DataFrame()
    df.loc[:, 'sample_id'] = source_df['SAMPID']
    df.loc[source_df.SAMPID, 'batch'] = 'GTEx_v8'
    df.loc[source_df.SAMPID, 'bam_path'] = source_df['rnaseq_bam']
    df.loc[source_df.SAMPID, 'bai_path'] = source_df['rnaseq_bai']
    df.loc[source_df.SAMPID, 'batch_detail'] = "GTEx_v8 " + source_df['SMNABTCH']
    #df.loc[source_df.SAMPID, 'output_dir'] = f"gs://macarthurlab-rnaseq/gtex_v8/fraser_count_rna/"
    df.loc[source_df.SAMPID, 'project'] = "gtex_v8"
    df.loc[source_df.SAMPID, 'sex'] = source_df['SEX']
    df.loc[source_df.SAMPID, 'age'] = source_df['AGE']
    df.loc[source_df.SAMPID, 'cause_of_death'] = source_df['DTHHRDY']
    df.loc[source_df.SAMPID, 'ancestry'] = None
    df.loc[source_df.SAMPID, 'tissue'] = source_df['SMTS']
    df.loc[source_df.SAMPID, 'tissue_detail'] = source_df['SMTSD']
    df.loc[source_df.SAMPID, 'read_length'] = source_df['SMRDLGTH']
    df.loc[source_df.SAMPID, 'stranded'] = "no"
    df.loc[source_df.SAMPID, 'RIN'] = source_df['SMRIN']

    return pd.concat([samples_df, df], axis="rows")


TISSUE_NAME_TO_SMTD = {
    "muscle": "Muscle - Skeletal",
    "fibroblasts": "Cells - Cultured fibroblasts",
    "whole_blood": "Whole Blood",
    "lymphocytes": "Cells - EBV-transformed lymphocytes",
}


def main():
    rnaseq_sample_metadata_df = get_joined_metadata_df()
    rnaseq_sample_metadata_df = rnaseq_sample_metadata_df[  # filter to samples
        rnaseq_sample_metadata_df.sample_id.isin(ANALYSIS_BATCHES["all_analysis_samples"]["samples"])
    ]
    gtex_rnaseq_sample_metadata_df = get_gtex_rnaseq_sample_metadata_df()
    #gtex_rnaseq_sample_metadata_df = gtex_rnaseq_sample_metadata_df #.rename({'SAMPID': 'sample_id'}, axis="columns").set_index('sample_id', drop=False)

    all_rdg_samples_df = None
    all_rdg_counts_df = None

    all_rdg_and_gtex_samples_df = None
    all_rdg_and_gtex_counts_df = None

    for tissue_name in ["whole_blood", "lymphocytes", "fibroblasts", "muscle"]:
        logging.info("----------------")
        logging.info(f"{tissue_name}:")
        samples_df = pd.DataFrame()
        SMTD_value = TISSUE_NAME_TO_SMTD[tissue_name]

        tgg_df = rnaseq_sample_metadata_df[rnaseq_sample_metadata_df['imputed tissue'] == tissue_name]
        logging.info(f"Got {len(tgg_df)} TGG samples for {tissue_name}")
        samples_df = transfer_metadata_columns_from_df(samples_df, tgg_df)
        samples_df['tissue'] = tissue_name
        samples_df['tissue_detail'] = SMTD_value
        if all_rdg_samples_df is None:
            all_rdg_samples_df = samples_df
        else:
            all_rdg_samples_df = pd.concat([all_rdg_samples_df, samples_df])

        gtex_df = gtex_rnaseq_sample_metadata_df[gtex_rnaseq_sample_metadata_df.SMTSD == SMTD_value]
        gtex_df = gtex_df.sort_values(by='SMRIN', ascending=False)
        logging.info(f"Got {len(gtex_df)} GTEx samples for {tissue_name}")

        if EXPORT_100_MALE_AND_100_FEMALE_GTEX_SAMPLES:
            samples_df = transfer_metadata_columns_from_GTEx_df(samples_df, gtex_df[gtex_df['SEX'] == 'M'][:100], tissue_name)
            samples_df = transfer_metadata_columns_from_GTEx_df(samples_df, gtex_df[gtex_df['SEX'] == 'F'][:100], tissue_name)
        else:
            samples_df = transfer_metadata_columns_from_GTEx_df(samples_df, gtex_df[:100], tissue_name)

        samples_df['tissue'] = tissue_name

        if all_rdg_and_gtex_samples_df is None:
            all_rdg_and_gtex_samples_df = samples_df
        else:
            all_rdg_and_gtex_samples_df = pd.concat([all_rdg_and_gtex_samples_df, samples_df])

        #tsv_output_path = f"metadata_table_for_{tissue_name}.tsv"
        #samples_df.to_csv(tsv_output_path, sep="\t", index=False)
        #print(f"Wrote {len(samples_df)} samples to {tsv_output_path}")

        print("TGG table columns:")
        print(tgg_df.columns)

        print("GTEx table columns:")
        print(gtex_df.columns)

        print("Output columns:")
        print(samples_df.columns)

        # create OUTRIDER data input table
        RDG_file_paths = {}
        for file_path in glob.glob(RDG_GENE_READS_PATH):
            sample_id = os.path.basename(file_path).replace(".gene_reads.gct.gz", "")
            if sample_id in set(samples_df.sample_id):
                RDG_file_paths[sample_id] = file_path

        RDG_sample_ids = [s for s in list(samples_df.sample_id) if "GTEX" not in s]
        GTEX_sample_ids = [s for s in list(samples_df.sample_id) if "GTEX" in s]

        if len(RDG_file_paths) != len(RDG_sample_ids):
            print(f"ERROR: len(file_paths) != len(samples_df): {len(RDG_file_paths)} != {len(samples_df)}" )
            continue
        print(f"Found all {len(RDG_sample_ids)} RDG samples and {len(GTEX_sample_ids)} GTEX samples for {tissue_name}")

        print(f"Reading {len(GTEX_sample_ids)} GTEX samples from {GTEX_GENE_READS_PATH}")
        with gzip.open(os.path.expanduser(GTEX_GENE_READS_PATH)) as f:
            next(f); next(f)
            gtex_counts_df = pd.read_table(f, usecols=['Name'] + list(GTEX_sample_ids))
        gtex_counts_df = gtex_counts_df.rename({'Name': 'gene_id'}, axis="columns")
        gtex_counts_df = gtex_counts_df.set_index('gene_id')

        rdg_counts_df = None
        for sample_id, file_path in RDG_file_paths.items():
            print(f"Reading {file_path}")
            with gzip.open(file_path, "rt") as f:
                next(f); next(f); next(f)
                current_gene_reads_df = pd.read_table(f, names=["gene_id", "name", sample_id], usecols=['gene_id', sample_id])
                current_gene_reads_df = current_gene_reads_df.set_index('gene_id')
                if rdg_counts_df is None:
                    rdg_counts_df = current_gene_reads_df
                else:
                    rdg_counts_df = pd.merge(rdg_counts_df, current_gene_reads_df, left_index=True, right_index=True, how="outer")

        print(len(rdg_counts_df.index), len(set(gtex_counts_df.index)), len(set(rdg_counts_df.index) & set(gtex_counts_df.index)))
        rdg_and_gtex_counts_df = pd.merge(rdg_counts_df, gtex_counts_df, left_index=True, right_index=True, how="outer")
        rdg_and_gtex_counts_df = rdg_and_gtex_counts_df.fillna(0)
        #tsv_output_path = f"OUTRIDER_input_table_for_{tissue_name}.tsv"
        #rdg_and_gtex_counts_df.reset_index().to_csv(tsv_output_path, sep="\t", index=False)
        #print(f"Wrote {len(rdg_and_gtex_counts_df)} genes x {len(rdg_and_gtex_counts_df.columns)} samples to {tsv_output_path}")

        if all_rdg_counts_df is None:
            all_rdg_counts_df = rdg_counts_df
        else:
            all_rdg_counts_df = pd.merge(all_rdg_counts_df, rdg_counts_df, left_index=True, right_index=True, how="outer")

        if all_rdg_and_gtex_counts_df is None:
            all_rdg_and_gtex_counts_df = rdg_and_gtex_counts_df
        else:
            all_rdg_and_gtex_counts_df = pd.merge(all_rdg_and_gtex_counts_df, rdg_and_gtex_counts_df, left_index=True, right_index=True, how="outer")

    #tsv_output_path = f"metadata_table_for_all_RDG_samples.tsv"
    #all_rdg_samples_df.to_csv(tsv_output_path, sep="\t", index=False)
    #print(f"Wrote {len(all_rdg_samples_df)} samples to {tsv_output_path}")

    tsv_output_path = f"metadata_table_for_all_RDG_and_GTEX_samples.tsv"
    all_rdg_and_gtex_samples_df.to_csv(tsv_output_path, sep="\t", index=False)
    print(f"Wrote {len(all_rdg_and_gtex_samples_df)} samples to {tsv_output_path}")
    #os.system(f"gsutil -m cp {tsv_output_path} gs://seqr-bw/project__rnaseq/")

    #all_rdg_counts_df = all_rdg_counts_df.fillna(0)
    #tsv_output_path = f"OUTRIDER_input_table_RDG_counts_for_all_tissues.tsv"
    #all_rdg_counts_df.reset_index().to_csv(tsv_output_path, sep="\t", index=False)
    #print(f"Wrote {len(all_rdg_counts_df)} genes x {len(all_rdg_counts_df.columns)} samples to {tsv_output_path}")

    all_rdg_and_gtex_counts_df = all_rdg_and_gtex_counts_df.fillna(0)
    tsv_output_path = f"OUTRIDER_input_table_RDG_and_GTEX_counts_for_all_tissues.tsv.gz"
    all_rdg_and_gtex_counts_df.reset_index().to_csv(tsv_output_path, sep="\t", index=False)
    print(f"Wrote {len(all_rdg_and_gtex_counts_df)} genes x {len(all_rdg_and_gtex_counts_df.columns)} samples to {tsv_output_path}")

    #os.system(f"gsutil -m cp {tsv_output_path} gs://seqr-bw/project__rnaseq/")


if __name__ == "__main__":
    main()
