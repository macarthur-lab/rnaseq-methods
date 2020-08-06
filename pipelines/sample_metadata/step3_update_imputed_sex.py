
#%%
import io
import math
from pprint import pprint
import pandas as pd
import hail as hl
import gzip
import requests
import traceback

from gspread_dataframe import set_with_dataframe

from sample_metadata.utils import \
    get_joined_metadata_df, \
    get_seqr_info_and_other_metadata_df, \
    get_seqr_info_and_other_metadata_worksheet


#%%
data_paths_df_with_metadata_columns = get_joined_metadata_df()
data_paths_df_with_metadata_columns.columns

#%%
sex_biased_genes_df = pd.read_table("https://raw.githubusercontent.com/berylc/MendelianRNA-seq/master/data/sex_biased_genes.txt")
sex_biased_genes_df['gene_id'] = sex_biased_genes_df['gene_id'].apply(lambda x: x.split(".")[0])
sex_biased_genes_df = sex_biased_genes_df.set_index('gene_id', drop=False)
sex_biased_genes_df = sex_biased_genes_df[['gene_id', 'gene_name', 'GeneGroup', 'chr', 'coeff']]

sex_biased_genes_df

#%%

imputed_sex_list = []
imputed_sex_doesnt_match_list = []
for _, row in data_paths_df_with_metadata_columns.iterrows():
    path = row['rnaseqc_gene_reads']
    sample_id = row['sample_id']

    print(sample_id, path)
    try:
        previous_sex = (row['sex'] if row['sex'] and not isinstance(row['sex'], float) else '').strip()
        imputed_sex = (row['imputed sex'] if row['imputed sex'] and not isinstance(row['imputed sex'], float) else '').strip()
        weighted_sum = row['weighted_sex_gene_expression']

        if not imputed_sex:
            file = hl.hadoop_open(row['rnaseqc_gene_reads'], "r")
            #if path.endswith(".gz"):
            #    file = gzip.GzipFile(fileobj=file)

            next(file) # skip header
            next(file)
            next(file)

            gene_reads_df = pd.read_table(file, names=["gene_id", "name", "count"])
            gene_reads_df['gene_id'] = gene_reads_df['gene_id'].apply(lambda x: x.split(".")[0])
            gene_reads_df = gene_reads_df.set_index('gene_id', drop=False)
            shared_gene_ids = set(sex_biased_genes_df.gene_id) & set(gene_reads_df.gene_id)
            weighted_sum = sum(gene_reads_df.loc[shared_gene_ids]['count'] * sex_biased_genes_df.loc[shared_gene_ids]['coeff'])
            imputed_sex = "" if abs(weighted_sum) < 500 else ("M" if weighted_sum < 0 else "F")

        values = {
            'sample_id': sample_id,
            'sex': previous_sex,
            'imputed sex': imputed_sex,
            'weighted_sex_gene_expression': weighted_sum,
        }
        print(", ".join(map(str, values.values())))

        if previous_sex and previous_sex != imputed_sex:
            print("WARNING: previous != imputed", previous_sex, imputed_sex)
            imputed_sex_doesnt_match_list.append(values)

        imputed_sex_list.append(values)

        file.close()

    except Exception as e:
        print(e)
        pprint(dict(row))
        traceback.print_exc()


#%%

df = get_seqr_info_and_other_metadata_df()
df = df.set_index('sample_id', drop=False)

df.columns

#%%

imputed_sex_df = pd.DataFrame(imputed_sex_list)
imputed_sex_df = imputed_sex_df.set_index('sample_id')

imputed_sex_df.columns

assert set(imputed_sex_df.index) == set(df.sample_id), set(imputed_sex_df.index) - set(df.sample_id)

df['imputed sex'] = imputed_sex_df['imputed sex']
df['weighted_sex_gene_expression'] = imputed_sex_df['weighted_sex_gene_expression']

df.columns

#%%
# export joined data to SEQR_INFO_AND_OTHER_METADATA_WORKSHEET
ws = get_seqr_info_and_other_metadata_worksheet()
set_with_dataframe(ws, df.fillna(''), resize=True)

print("Updated", ws.title)


#%%
print("---------------")
for v in imputed_sex_doesnt_match_list:
    print(v)
