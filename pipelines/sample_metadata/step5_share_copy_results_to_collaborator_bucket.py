import datetime
import logging
import pandas as pd

from sample_metadata.rnaseq_metadata_utils import get_rnaseq_metadata_joined_with_paths_df, \
    get_gtex_rnaseq_sample_metadata_df, \
    get_rnaseq_downstream_analysis_metadata_worksheet

from gspread_dataframe import set_with_dataframe

logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s', level=logging.INFO)
logger = logging.getLogger(__name__)

#%%

metadata_df = get_rnaseq_metadata_joined_with_paths_df()

#%%

# filter
metadata_df = metadata_df[metadata_df.sample_id.str.startswith("BON_")]

destination_bucket = "gs://tgg-rnaseq-bonnemann/"
#%%

metadata_df.columns

import hail as hl
hl.init(log="/dev/null")

#hl.hadoop_copy("gs://tgg-rnaseq-bonnemann/")


#%%

metadata_df = metadata_df.drop(columns=['sex']).rename(columns={
    'star_pipeline_batch': 'batch',
    'batch_date_from_hg19_bam_header': 'sequencing_date',
    'star_bam': 'bam_path',
    'star_bai': 'bai_path',
    'imputed sex': 'sex',
    'imputed tissue': 'tissue',
    'read length (rnaseqc)': 'read_length',
    'stranded? (rnaseqc)': 'stranded',
    'coverage_bigwig': 'bigWig',
    'grch38_vcf': 'vcf',
})

wgs_projects = metadata_df['proj WGS (seqr)'].apply(lambda p: (p and p.strip()) or None)
wes_projects = metadata_df['proj WES (seqr)'].apply(lambda p: p or '')
metadata_df.loc[:, 'project'] = wgs_projects.fillna(wes_projects)

#%%
metadata_df.loc[:, 'age'] = metadata_df['Age at muscle biopsy (Beryl:Supp.)'].apply(lambda x: ('' if not x or not x.strip() else ('0' if (x and (x == 'At birth' or 'm' in x)) else x.replace('y', ''))))
metadata_df.loc[:, 'tissue_detail'] = metadata_df.tissue.apply(lambda t: TISSUE_NAME_TO_SMTSD[t])
metadata_df.loc[:, 'is_GTEx_sample'] = False

for empty_column in 'cause_of_death', 'RIN':  #'ancestry',
    metadata_df.loc[:, empty_column] = None

#%%
metadata_df = metadata_df[OUTPUT_COLUMNS]

#%%

## add GTEx
gtex_df = get_gtex_rnaseq_sample_metadata_df()
gtex_df = gtex_df[gtex_df.SMTSD.isin(TISSUE_NAME_TO_SMTSD.values())]

print(f"Got {len(gtex_df)} GTEx samples")
#%%

gtex_df = gtex_df.rename(columns={
    'SAMPID': 'sample_id',
    'rnaseq_bam': 'bam_path',
    'rnaseq_bai': 'bai_path',
    'SEX': 'sex',
    'AGE': 'age',
    'DTHHRDY': 'cause_of_death',
    'SMTSD': 'tissue_detail',
    'SMRDLGTH': 'read_length',
    'SMRIN': 'RIN',
    'wgs_vcf': 'vcf',
})

#%%

#sequencing_date = get_date_from_bam_header(next(iter(gtex_df.bam_path)))
sequencing_date = '2014-02'  #  this is the date from one of the GTEx bams

#%%

gtex_df.loc[:, 'tissue'] = gtex_df.tissue_detail.apply(lambda t: SMTSD_TO_TISSUE_NAME[t])
gtex_df.loc[:, 'batch'] = 'GTEx_v8'
gtex_df.loc[:, 'project'] = 'gtex_v8'
gtex_df.loc[:, 'stranded'] = 'no'
gtex_df.loc[:, 'sequencing_date'] = sequencing_date
gtex_df.loc[:, 'project'] = "gtex_v8"
#gtex_df.loc[:, 'ancestry'] = None
gtex_df.loc[:, 'is_GTEx_sample'] = True

#'star_SJ_out_tab', 'bigWig', 'wgs_vcf', 'wes_vcf'

#%%
gtex_df = gtex_df[OUTPUT_COLUMNS]

#%%

for column in gtex_df.columns:
    print('-'*100)
    print(column)
    print(gtex_df[column])

#%%

for tissue_smtsd in TISSUE_NAME_TO_SMTSD.values():
    gtex_df_for_tissue = gtex_df[gtex_df.tissue_detail == tissue_smtsd]
    gtex_df_for_tissue.sort_values(by='RIN', ascending=False, inplace=True)
    print(f"Found {len(gtex_df_for_tissue)} GTEx {tissue_smtsd} samples")
    metadata_df = pd.concat([metadata_df, gtex_df_for_tissue[:100]])

#df.loc[source_df.sample_id, 'output_dir'] = source_df['star_pipeline_batch'].apply(lambda batch_name: f"gs://tgg-rnaseq/{batch_name}/fraser_count_rna/")
#df.loc[source_df.SAMPID, 'output_dir'] = f"gs://tgg-rnaseq/gtex_v8/fraser_count_rna/"

#%%
metadata_df.sort_values(by=['is_GTEx_sample', 'batch', 'tissue', 'project'], ascending=True, inplace=True)

ws = get_rnaseq_downstream_analysis_metadata_worksheet()
set_with_dataframe(ws, metadata_df.fillna(''), resize=True)

print("Updated", ws.title)

#%%

timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
tsv_output_path = f"~/project__rnaseq/code/rnaseq_methods/pipelines/sample_metadata/rnaseq_downstream_analysis_metadata__{timestamp}.tsv"
metadata_df.to_csv(tsv_output_path, sep="\t", index=False)
print(f"Wrote {len(metadata_df)} samples to {tsv_output_path}")


#%%
print("Sample counts by batch:")
metadata_df.groupby(["batch", "tissue"])["sample_id"].nunique()

#%%
print("GTEx sample counts:")
metadata_df[metadata_df["batch"] == "GTEx_v8"].groupby(["tissue"])["sample_id"].nunique()

#%%
print("Non-GTEx sample counts:")
metadata_df[metadata_df["batch"] != "GTEx_v8"].groupby(["tissue"])["sample_id"].nunique()

#%%
