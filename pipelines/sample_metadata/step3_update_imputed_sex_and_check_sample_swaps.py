
#%%
import datetime
import hail as hl
import os
from pprint import pprint
import pandas as pd
import traceback

from gspread_dataframe import set_with_dataframe

from sample_metadata.rnaseq_metadata_utils import \
    get_rnaseq_metadata_joined_with_paths_df, \
    get_rnaseq_metadata_df, \
    get_rnaseq_metadata_worksheet

#%%
hl.init(log="/dev/null")

#%%
rnaseq_metadata_joined_with_paths_df = get_rnaseq_metadata_joined_with_paths_df()
rnaseq_metadata_joined_with_paths_df.columns


#%%

# Check for sample swaps:
#
os.system(
    "gsutil cat gs://tgg-rnaseq/*/somalier_sample_swap/*.groups.tsv " +
    "| awk '{ if($2 < 0.75) { print( $2, $1 ) } }' " +
    "| sed 's/,sm1//g' " +
    "| sort -n -r")

# TODO rewrite this in python

#%%
#sex_biased_genes_df = pd.read_table("https://raw.githubusercontent.com/berylc/MendelianRNA-seq/master/data/sex_biased_genes.txt")
#sex_biased_genes_df['gene_id'] = sex_biased_genes_df['gene_id'].apply(lambda x: x.split(".")[0])
#sex_biased_genes_df = sex_biased_genes_df.set_index('gene_id', drop=False)
#sex_biased_genes_df = sex_biased_genes_df[['gene_id', 'gene_name', 'GeneGroup', 'chr', 'coeff']]

#sex_biased_genes_df

#%%

imputed_sex_list = []
imputed_sex_doesnt_match_list = []

for _, row in rnaseq_metadata_joined_with_paths_df.iterrows():
    sample_id = row['sample_id']

    try:
        previous_sex = (row['sex'] if row['sex'] and not isinstance(row['sex'], float) else '').strip()
        imputed_sex = (row['imputed sex'] if row['imputed sex'] and not isinstance(row['imputed sex'], float) else '').strip()
        if True or not imputed_sex:
            somalier_results_path = f"gs://tgg-rnaseq/{row['star_pipeline_batch']}/somalier_sample_swap/{row['sample_id']}.somalier_results.tsv"
            somalier_results_df = pd.read_table(hl.hadoop_open(somalier_results_path, "r"))
            if "chrX/chrY" not in set(somalier_results_df.columns):
                print(f"ERROR: chrX/chrY not found in {somalier_results_path}")
                imputed_sex = ""
            else:
                coverage_ratio = somalier_results_df["chrX/chrY"].loc[0]
                imputed_sex = "M" if coverage_ratio < 100 else "F"

            #file = hl.hadoop_open(row['rnaseqc_gene_reads'], "r")
            #if path.endswith(".gz"):
            #    file = gzip.GzipFile(fileobj=file)

            #next(file) # skip header
            #next(file)
            #next(file)

            #gene_reads_df = pd.read_table(file, names=["gene_id", "name", "count"])
            #file.close()

            #gene_reads_df['gene_id'] = gene_reads_df['gene_id'].apply(lambda x: x.split(".")[0])
            #gene_reads_df = gene_reads_df.set_index('gene_id', drop=False)
            #shared_gene_ids = set(sex_biased_genes_df.gene_id) & set(gene_reads_df.gene_id)
            #weighted_sum = sum(gene_reads_df.loc[shared_gene_ids]['count'] * sex_biased_genes_df.loc[shared_gene_ids]['coeff'])
            #weighted_sum *= 10.0**5/sum(gene_reads_df['count'])

        values = {
            'sample_id': sample_id,
            'sex': previous_sex,
            'imputed sex': imputed_sex,
            'chrX/chrY coverage': coverage_ratio,
            #'weighted_sex_gene_expression': weighted_sum,
        }

        print(", ".join(map(str, values.values())))

        if previous_sex and previous_sex != imputed_sex:
            print("WARNING: previous != imputed", previous_sex, imputed_sex)
            imputed_sex_doesnt_match_list.append(values)

        imputed_sex_list.append(values)

    except Exception as e:
        print(e)
        pprint(dict(row))
        traceback.print_exc()


#%%

df = get_rnaseq_metadata_df()
df = df.set_index('sample_id', drop=False)

df.columns

#%%

imputed_sex_df = pd.DataFrame(imputed_sex_list)
imputed_sex_df = imputed_sex_df.set_index('sample_id')

imputed_sex_df.columns

assert set(imputed_sex_df.index) == set(df.sample_id), set(imputed_sex_df.index) - set(df.sample_id)

df['imputed sex'] = imputed_sex_df['imputed sex']
#df['chrX/chrY coverage'] = imputed_sex_df['chrX/chrY coverage']
#df['weighted_sex_gene_expression'] = imputed_sex_df['weighted_sex_gene_expression']

df.columns


#%%
print("---------------")
for v in imputed_sex_doesnt_match_list:
    print(v)

#%%
# export joined data to RNASEQ_METADATA_WORKSHEET
ws = get_rnaseq_metadata_worksheet()
set_with_dataframe(ws, df.fillna(''), resize=True)

print("Updated", ws.title)


#%%

timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
tsv_output_path = f"~/project__rnaseq/code/rnaseq_methods/pipelines/sample_metadata/rnaseq_metadata_with_sex_and_sample_swaps__{timestamp}.tsv"
df.to_csv(tsv_output_path, sep="\t", index=False)
print(f"Wrote {len(df)} samples to {tsv_output_path}")

#%%


