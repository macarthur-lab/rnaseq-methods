import datetime
import logging
import pandas as pd

from sample_metadata.rnaseq_metadata_utils import get_rnaseq_downstream_analysis_metadata_df

logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s', level=logging.INFO)
logger = logging.getLogger(__name__)

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

OUTPUT_COLUMNS = ['DROP_GROUP', 'RNA_ID', 'DNA_ID', 'RNA_BAM_FILE', 'RNA_BAI_FILE', 'DNA_VCF_FILE', 'PAIRED_END', 'STRAND', 'HPO_TERMS', 'COUNT_MODE', 'COUNT_OVERLAPS', 'ANNOTATION', 'GENE_COUNTS_FILE']

#%%
# get metadata table and filter out excluded samples
rnaseq_sample_metadata_df = get_rnaseq_downstream_analysis_metadata_df()
all_rdg_and_gtex_samples_df = None
for tissue_name in ["whole_blood", "lymphocytes", "fibroblasts", "muscle"]:
    logging.info("----------------")
    logging.info(f"{tissue_name}:")
    samples_df = rnaseq_sample_metadata_df[rnaseq_sample_metadata_df.tissue == tissue_name]

    # process RDG samples
    samples_df_rdg = samples_df[~samples_df.sample_id.str.startswith('GTEX')]
    samples_df_rdg.loc[:, 'DROP_GROUP'] = f"{tissue_name},{tissue_name}_with_GTEx" if len(samples_df_rdg) > 30 else f"{tissue_name}_with_GTEx"
    samples_df_rdg.loc[:, 'RNA_ID'] = samples_df_rdg['sample_id']
    samples_df_rdg.loc[:, 'DNA_ID'] = samples_df_rdg['sample_id']
    samples_df_rdg.loc[:, 'RNA_BAM_FILE'] = samples_df_rdg['bam_path']
    samples_df_rdg.loc[:,  'RNA_BAI_FILE'] = samples_df_rdg['bai_path']
    samples_df_rdg.loc[:,  'DNA_VCF_FILE'] = samples_df_rdg['vcf']
    samples_df_rdg.loc[:,  'PAIRED_END'] = True
    samples_df_rdg.loc[:,  'STRAND'] = samples_df_rdg['stranded']  #  "yes", "no"
    samples_df_rdg.loc[:,  'HPO_TERMS'] = ""
    samples_df_rdg.loc[:,  'COUNT_MODE'] = "IntersectionStrict"
    samples_df_rdg.loc[:,  'COUNT_OVERLAPS'] = True
    samples_df_rdg.loc[:,  'ANNOTATION'] = "v34"
    samples_df_rdg.loc[:,  'GENE_COUNTS_FILE'] = ""
    samples_df_rdg = samples_df_rdg[OUTPUT_COLUMNS]

    logging.info(f"Got {len(samples_df_rdg)} TGG samples for {tissue_name}")

    # process GTEx samples
    samples_df_gtex = samples_df[samples_df.is_GTEx_sample == "TRUE"]
    #samples_df_gtex = samples_df_gtex[samples_df_gtex["vcf"].str.len() > 1]
    samples_df_gtex.loc[:, 'DROP_GROUP'] = f"{tissue_name}_with_GTEx"
    samples_df_gtex.loc[:, 'RNA_ID'] = samples_df_gtex['sample_id']
    samples_df_gtex.loc[:, 'DNA_ID'] = samples_df_gtex['sample_id']
    samples_df_gtex.loc[:, 'RNA_BAM_FILE'] = samples_df_gtex['bam_path']
    samples_df_gtex.loc[:, 'RNA_BAI_FILE'] = samples_df_gtex['bai_path']
    #samples_df_gtex.loc[:, 'DNA_VCF_FILE'] = samples_df_gtex['wes_vcf']
    samples_df_gtex.loc[:, 'DNA_VCF_FILE'] = samples_df_gtex['vcf']
    samples_df_gtex.loc[:, 'PAIRED_END'] = True
    samples_df_gtex.loc[:, 'STRAND'] = "no"
    samples_df_gtex.loc[:, 'HPO_TERMS'] = ""
    samples_df_gtex.loc[:, 'COUNT_MODE'] = "IntersectionStrict"
    samples_df_gtex.loc[:, 'COUNT_OVERLAPS'] = True
    samples_df_gtex.loc[:, 'ANNOTATION'] = "v34"
    samples_df_gtex.loc[:, 'GENE_COUNTS_FILE'] = ""

    logging.info(f"Got {len(samples_df_gtex)} GTEx samples for {tissue_name}")

    if all_rdg_and_gtex_samples_df is None:
        all_rdg_and_gtex_samples_df = pd.concat([samples_df_rdg, samples_df_gtex])
    else:
        all_rdg_and_gtex_samples_df = pd.concat([all_rdg_and_gtex_samples_df, samples_df_rdg, samples_df_gtex])

    #tsv_output_path = f"metadata_table_for_{tissue_name}.tsv"
    #samples_df.to_csv(tsv_output_path, sep="\t", index=False)
    #print(f"Wrote {len(samples_df)} samples to {tsv_output_path}")

# create separate groups for MAE analysis, including only samples that have a VCF
#all_rdg_and_gtex_samples_df.loc[
#    all_rdg_and_gtex_samples_df["DNA_VCF_FILE"].str.len() > 1,
#    'DROP_GROUP'] = all_rdg_and_gtex_samples_df['DROP_GROUP'].apply(
#        lambda groups: groups+","+",".join([f"{g}_MAE" for g in groups.split(",")]))

#%%
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

#%%
