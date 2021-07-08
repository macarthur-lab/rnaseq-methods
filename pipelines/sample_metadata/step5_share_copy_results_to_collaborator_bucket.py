import logging
import os
import pandas as pd


from sample_metadata.rnaseq_metadata_utils import get_rnaseq_metadata_joined_with_paths_df, \
    get_rnaseq_downstream_analysis_metadata_worksheet

from gspread_dataframe import set_with_dataframe

logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s', level=logging.INFO)
logger = logging.getLogger(__name__)

metadata_df = get_rnaseq_metadata_joined_with_paths_df()
metadata_df = metadata_df[metadata_df.star_pipeline_batch == "batch_1_muntoni"]   # filter
list(metadata_df.columns)


metadata_df = metadata_df[[
    'sample_id',
    'star_pipeline_batch',
    'star_bam',
    'star_bai',
    'star_SJ_out_tab',
    'star_reads_per_gene_tab',
    'grch38_vcf',
    'rnaseqc_gene_reads',
    'rnaseqc_exon_reads',
    'rnaseqc_gene_tpm',
    'rnaseqc_metrics',
    'junctions_bed',
    'coverage_bigwig',
    'portcullis_filtered',
    'portcullis_all',
    'hg19_bam',
    'hg19_bai',
    'fastqc_zip',
    'batch_date_from_hg19_bam_header',
    'imputed tissue',
    'imputed sex',
     'stranded? (rnaseqc)',
    'read length (rnaseqc)',
    'total reads x 10^6 (rnaseqc)',
    'mapping rate (rnaseqc)',
    'indiv (seqr)',
    'proj WGS (seqr)',
    'fam WGS (seqr)',
    'proj WES (seqr)',
    'fam WES (seqr)',
    'genome (seqr)',
    'population (seqr)',
    'sample id (seqr)',
    'sample type (seqr)',
    'analysis status (seqr)',
    'variant tags (seqr)',
    'variant notes (seqr)',
    'coded phenotype (seqr)',
    'analysis summary + notes (seqr)',
]]

output_path = "/tmp/muntoni_rnaseq_sample_metadata.tsv"
metadata_df.to_csv(output_path, header=True, index=False, sep="\t")
print(f"Wrote {len(metadata_df)} rows to {output_path}")


#%%
#user = "svetlana.gorokhova@univ-amu.fr"
user = "veronica.pini01@ateneopv.it"
rows = [row for _, row in metadata_df.iterrows()]
batch_size = 20
for i in range(0, len(rows), batch_size):
    paths = []
    for row in rows[i:i+batch_size]:
        for key, value in row.to_dict().items():
            if value and value.startswith("gs://"):
                value = os.path.splitext(value)[0] + ".*"  # remove extension, add wildcard to also capture .tbi, etc.
                paths.append(value)
    os.system(f"gsutil acl ch -u {user}:R " + " ".join(paths))

# gsutil acl ch -u allUsers:R gs://tgg-viewer-configs/*.json

#%%


#import hail as hl
#hl.init(log="/dev/null")

#hl.hadoop_copy("gs://tgg-rnaseq-bonnemann/")

