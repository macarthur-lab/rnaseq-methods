
#%%
import argparse
import datetime
import glob
import gzip
import hashlib
import logging
import os
import pandas as pd
import sys

from sample_metadata.utils import get_joined_metadata_df, get_gtex_rnaseq_sample_metadata_df

logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s', level=logging.INFO)
logger = logging.getLogger(__name__)

#%%

df = get_joined_metadata_df()
#%%

df.columns

#%%

df[['proj (seqr)', 'imputed tissue']].reset_index().groupby(('imputed tissue')).count()

#%%

df = df[['proj (seqr)', 'imputed tissue', 'imputed sex', 'read length (rnaseqc)', 'batch_date_from_hg19_bam_header']].reset_index()

#%%

df = df[(df['imputed tissue'].str.len() > 2) & (df['batch_date_from_hg19_bam_header'].str.len() > 2)]
#%%

#df['year'] = df['batch_date_from_hg19_bam_header'].apply(lambda x: None if not x.strip() else int(x.split("-")[0]) > 2017)

#%%

df.groupby(['imputed tissue', 'imputed sex', 'read length (rnaseqc)']).count()

#%%

df.groupby('imputed tissue').count()

#%%
