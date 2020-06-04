
#%%

from __future__ import print_function
import os
import pandas as pd
import sys

from gspread_dataframe import set_with_dataframe

from sample_metadata.utils import get_seqr_info_and_other_metadata_df, SEQR_INFO_AND_OTHER_METADATA_WORKSHEET


#%%
df = get_seqr_info_and_other_metadata_df()

df

#%%

imputed_tissues_df = pd.read_table(os.path.expanduser("~/project__rnaseq/data/samples/expression/rnaseqc_tpm/imputed_tissues.tsv"))

imputed_tissues_df

#%%

assert set(imputed_tissues_df.SAMPID) == set(df.sample_id)

#%%

df.loc[:, 'imputed tissue'] = imputed_tissues_df.imputed_tissue

df['imputed tissue']

#%%

df[['sample_id', 'imputed tissue']]

#%%

df
#%%
# export joined data to SEQR_INFO_AND_OTHER_METADATA_WORKSHEET
set_with_dataframe(SEQR_INFO_AND_OTHER_METADATA_WORKSHEET, df.fillna(''), resize=True)

print("Updated", SEQR_INFO_AND_OTHER_METADATA_WORKSHEET.title)


# %%
