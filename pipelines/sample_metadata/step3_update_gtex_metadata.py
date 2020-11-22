
#%%
import os
import pandas as pd
import subprocess
from gspread_dataframe import set_with_dataframe
from sample_metadata.rnaseq_metadata_utils import \
    get_gtex_rnaseq_sample_metadata_worksheet, \
    get_gtex_wgs_sample_metadata_worksheet, \
    get_gtex_wes_sample_metadata_worksheet, \
    get_gtex_individual_metadata_worksheet
from google.cloud import storage

#%%

pd.set_option('display.max_columns', 500)
#%%

# sample info table
df_samples = pd.read_table(
    "https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt")

# add column for individual ids
df_samples['SUBJID'] = df_samples.SAMPID.apply(lambda s: "-".join(s.split("-")[0:2]))

#%%

# individual info table - has 4 columns: SUBJID   SEX    AGE  DTHHRDY
df_indivs = pd.read_table(
    "https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt")

df_indivs['SEX'] = df_indivs.SEX.apply(lambda s: "M" if s == 1 else ('F' if s == 2 else 'U'))

assert set(df_indivs.SEX == {'F', 'M'})

ws = get_gtex_individual_metadata_worksheet()
set_with_dataframe(ws, df_indivs.fillna(''), resize=True)

#%%

df_rnaseq_samples = df_samples[df_samples.SMAFRZE == "RNASEQ"].merge(df_indivs, how="left", on="SUBJID")
df_rnaseq_samples = df_rnaseq_samples.set_index('SAMPID', drop=False).dropna(axis='columns', how="all").fillna('')

df_wes_samples = df_samples[df_samples.SMAFRZE == "WES"].merge(df_indivs, how="left", on="SUBJID")
df_wes_samples = df_wes_samples.set_index('SAMPID', drop=False).dropna(axis='columns', how="all").fillna('')

df_wgs_samples = df_samples[df_samples.SMAFRZE == "WGS"].merge(df_indivs, how="left", on="SUBJID").dropna(axis='columns', how="all").fillna('')
df_wgs_samples = df_wgs_samples.set_index('SAMPID', drop=False).dropna(axis='columns', how="all").fillna('')


#%%

# https://app.terra.bio/#workspaces/anvil-datastorage/AnVIL_GTEx_V8_hg38

"""
GTEx bucket contains:   $ gs -u seqr-project ls gs://fc-secure-ff8156a3-ddf3-42e4-9211-0fd89da62108/
gs://fc-secure-ff8156a3-ddf3-42e4-9211-0fd89da62108/GTEx_Analysis_2017-06-05_v8_ASE_WASP_chrX_raw_counts_by_subject/
gs://fc-secure-ff8156a3-ddf3-42e4-9211-0fd89da62108/GTEx_Analysis_2017-06-05_v8_ASE_WASP_counts_by_subject/
gs://fc-secure-ff8156a3-ddf3-42e4-9211-0fd89da62108/GTEx_Analysis_2017-06-05_v8_ASE_chrX_raw_counts_by_subject/
gs://fc-secure-ff8156a3-ddf3-42e4-9211-0fd89da62108/GTEx_Analysis_2017-06-05_v8_ASE_counts_by_subject/
gs://fc-secure-ff8156a3-ddf3-42e4-9211-0fd89da62108/GTEx_Analysis_2017-06-05_v8_Annotations/
gs://fc-secure-ff8156a3-ddf3-42e4-9211-0fd89da62108/GTEx_Analysis_2017-06-05_v8_ChIPseq_BAM_files/
gs://fc-secure-ff8156a3-ddf3-42e4-9211-0fd89da62108/GTEx_Analysis_2017-06-05_v8_DroNcseq_FASTQ_files/
gs://fc-secure-ff8156a3-ddf3-42e4-9211-0fd89da62108/GTEx_Analysis_2017-06-05_v8_Histology_SVS_files/
gs://fc-secure-ff8156a3-ddf3-42e4-9211-0fd89da62108/GTEx_Analysis_2017-06-05_v8_RIPseq_FASTQ_files/
gs://fc-secure-ff8156a3-ddf3-42e4-9211-0fd89da62108/GTEx_Analysis_2017-06-05_v8_RNA_MuTect/
gs://fc-secure-ff8156a3-ddf3-42e4-9211-0fd89da62108/GTEx_Analysis_2017-06-05_v8_RNAseq_BAM_files/
gs://fc-secure-ff8156a3-ddf3-42e4-9211-0fd89da62108/GTEx_Analysis_2017-06-05_v8_RNAseq_LeafCutter_junc_files/
gs://fc-secure-ff8156a3-ddf3-42e4-9211-0fd89da62108/GTEx_Analysis_2017-06-05_v8_RNAseq_aux_files/
gs://fc-secure-ff8156a3-ddf3-42e4-9211-0fd89da62108/GTEx_Analysis_2017-06-05_v8_RNAseq_bigWig_files/
gs://fc-secure-ff8156a3-ddf3-42e4-9211-0fd89da62108/GTEx_Analysis_2017-06-05_v8_WES_BAM_files/
gs://fc-secure-ff8156a3-ddf3-42e4-9211-0fd89da62108/GTEx_Analysis_2017-06-05_v8_WES_VCF_files/
gs://fc-secure-ff8156a3-ddf3-42e4-9211-0fd89da62108/GTEx_Analysis_2017-06-05_v8_WGS_CRAM_files/
gs://fc-secure-ff8156a3-ddf3-42e4-9211-0fd89da62108/GTEx_Analysis_2017-06-05_v8_WGS_VCF_files/
gs://fc-secure-ff8156a3-ddf3-42e4-9211-0fd89da62108/GTEx_Analysis_2017-06-05_v8_phASER/
gs://fc-secure-ff8156a3-ddf3-42e4-9211-0fd89da62108/GTEx_Analysis_2017-06-05_v8_phASER_WASP/
"""


def get_gtex_file_paths(directory):
    bucket = storage.Client().get_bucket('fc-secure-ff8156a3-ddf3-42e4-9211-0fd89da62108')
    gtex_bucket_blobs = bucket.list_blobs(prefix=directory)
    for b in gtex_bucket_blobs:
        yield b.public_url.replace("https://storage.googleapis.com/", "gs://")


def get_sample_id_from_path(p):
    return os.path.basename(p).split(".")[0]

#%%

# RNA-seq hg38 .bam, .bam.bai files
gtex_hg38_rnaseq_bams = {get_sample_id_from_path(p): p for p in get_gtex_file_paths("GTEx_Analysis_2017-06-05_v8_RNAseq_BAM_files") if p.endswith(".bam")}
print("Found %s hg38 RNA-seq bams" % len(gtex_hg38_rnaseq_bams))

df_rnaseq_samples.loc[:, 'rnaseq_bam'] = pd.Series(gtex_hg38_rnaseq_bams)
df_rnaseq_samples.loc[:, 'rnaseq_bai'] = pd.Series({k: v.replace('.bam', '.bam.bai') for k, v in gtex_hg38_rnaseq_bams.items()})

#%%

gtex_star_SJ_out_tab = {get_sample_id_from_path(p): p for p in get_gtex_file_paths("GTEx_Analysis_2017-06-05_v8_RNAseq_aux_files") if p.endswith(".SJ.out.tab")}
print("Found %s star SJ.out.tab files" % len(gtex_star_SJ_out_tab))

df_rnaseq_samples.loc[:, 'star_SJ_out_tab'] = pd.Series(gtex_star_SJ_out_tab)

#%%

gtex_bigWig = {get_sample_id_from_path(p): p for p in get_gtex_file_paths("GTEx_Analysis_2017-06-05_v8_RNAseq_bigWig_files") if p.endswith(".bigWig")}
print("Found %s bigWig files" % len(gtex_bigWig))

df_rnaseq_samples.loc[:, 'bigWig'] = pd.Series(gtex_bigWig)

#%%

# WGS hg38 .cram, .crai files
gtex_hg38_WGS_crams = {get_sample_id_from_path(p): p for p in get_gtex_file_paths("GTEx_Analysis_2017-06-05_v8_WGS_CRAM_files") if p.endswith(".cram")}
print("Found %s hg38 WGS crams" % len(gtex_hg38_WGS_crams))

df_wgs_samples.loc[:, 'wgs_cram'] = pd.Series(gtex_hg38_WGS_crams)
df_wgs_samples.loc[:, 'wgs_crai'] = pd.Series({k: v.replace('.cram', '.crai') for k, v in gtex_hg38_rnaseq_bams.items()})

#gtex_hg38_WGS_vcf = "gs://fc-secure-ff8156a3-ddf3-42e4-9211-0fd89da62108/GTEx_Analysis_2017-06-05_v8_WGS_VCF_files/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.vcf.gz"

#%%

df_rnaseq_samples = df_rnaseq_samples.set_index("SUBJID", drop=False)

# WGS vcf
WGS_VCF_PATH = "gs://fc-secure-ff8156a3-ddf3-42e4-9211-0fd89da62108/GTEx_Analysis_2017-06-05_v8_WGS_VCF_files/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_866Indiv.vcf.gz"
subj_ids = subprocess.check_output(f"gsutil -u seqr-project cat {WGS_VCF_PATH} | gunzip -c - | head -n 1000 | grep CHROM | cut -f 10- | tr '\t' '\n'", shell=True, encoding="UTF-8").strip().split()
df_rnaseq_samples.loc[:, 'wgs_vcf'] = pd.Series({subj_id: WGS_VCF_PATH for subj_id in subj_ids})

# WES vcf
WES_VCF_PATH = "gs://fc-secure-ff8156a3-ddf3-42e4-9211-0fd89da62108/GTEx_Analysis_2017-06-05_v8_WES_VCF_files/GTEx_Analysis_2017-06-05_v8_WholeExomeSeq_979Indiv_VEP_annot.vcf.gz"
subj_ids = subprocess.check_output(f"gsutil -u seqr-project cat {WES_VCF_PATH} | gunzip -c - | head -n 1000 | grep CHROM | cut -f 10- | tr '\t' '\n'", shell=True, encoding="UTF-8").strip().split()
df_rnaseq_samples.loc[:, 'wes_vcf'] = pd.Series({subj_id: WES_VCF_PATH for subj_id in subj_ids})

#%%

ws = get_gtex_rnaseq_sample_metadata_worksheet()
set_with_dataframe(ws, df_rnaseq_samples, resize=True)

#%%

ws = get_gtex_wes_sample_metadata_worksheet()
set_with_dataframe(ws, df_wes_samples, resize=True)

#%%

ws = get_gtex_wgs_sample_metadata_worksheet()
set_with_dataframe(ws, df_wgs_samples, resize=True)


#%%
