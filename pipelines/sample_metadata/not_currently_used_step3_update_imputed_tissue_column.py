
#%%

from __future__ import print_function
import os
import pandas as pd
import sys

from gspread_dataframe import set_with_dataframe

from sample_metadata.rnaseq_metadata_utils import get_rnaseq_metadata_df, get_rnaseq_metadata_worksheet


#%%
df = get_rnaseq_metadata_df()
df = df.set_index('sample_id', drop=False)

df.columns
#%%

imputed_tissues_df = pd.read_table("~/project__rnaseq/data/samples/expression/rnaseqc_tpm/imputed_tissues.tsv")
imputed_tissues_df = imputed_tissues_df.set_index('SAMPID')
imputed_tissues_df = imputed_tissues_df.rename({"imputed_tissue": "imputed tissue"}, axis="columns")

imputed_tissues_df.columns

#%%

assert set(imputed_tissues_df.index) == set(df.sample_id)

#%%

df['imputed tissue'] = imputed_tissues_df['imputed tissue']


#%%
# export joined data to RNASEQ_METADATA_WORKSHEET
ws = get_rnaseq_metadata_worksheet()
set_with_dataframe(ws, df.fillna(''), resize=True)

print("Updated", ws.title)


# %%
