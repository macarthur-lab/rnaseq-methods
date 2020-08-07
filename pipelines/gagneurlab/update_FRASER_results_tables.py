#%%

import glob
import os
import pandas as pd
from pprint import pprint
from collections import defaultdict

#%%

# gsutil -m cp gs://macarthurlab-rnaseq/gagneur/outrider/results/*with*.tsv.gz .

all_tables = defaultdict(list)
results_tables = "/Users/weisburd/project__rnaseq/code/rnaseq_methods/pipelines/gagneurlab/FRASER_results1/using_AE/*.tsv.gz"
for path in [p for p in glob.glob(results_tables)]:
    name = os.path.basename(path)
    label = name.split("_with")[0]
    if "_without_GTEX" in name:
        has_GTEX = False
    elif "_with_GTEX" in name:
        has_GTEX = True
    else:
        raise ValueError(f"Unexpected filename: {name}")
    all_tables[label].append({'has_GTEX': has_GTEX, 'path': path})

pprint(sorted(all_tables))

#%%


def read_table(path, padj_threshold=0.10):
    df = pd.read_table(path)
    df = df[df.padjust < padj_threshold]
    df['geneID'] = df['geneID'].apply(lambda s: s.split(".")[0]) # remove gene versions
    df = df.set_index(["sampleID", "geneID"])
    return df


results = {}
for label, tables in all_tables.items():
    if len(tables) > 2:
        raise ValueError(f"More than 2 tables found for: {label}")

    if len(tables) == 1:
        results[label] = read_table(tables[0]['path']).reset_index()
        continue
    tables.sort(key=lambda x: x['has_GTEX'])
    t1 = read_table(tables[0]['path'])
    t2 = read_table(tables[1]['path'])

    #print(t1.columns)  # 'sampleID', 'geneID', 'pValue', 'padjust', 'zScore', 'rawcounts', 'q'
    result = t1.join(t2,
        how="outer",
        lsuffix=("_with_GTEX" if tables[0]['has_GTEX'] else "_without_GTEX"),
        rsuffix=("_with_GTEX" if tables[1]['has_GTEX'] else "_without_GTEX")).reset_index()
    results[label] = result


#%%

# https://docs.google.com/spreadsheets/d/1mvUcRDyAINf1utpF4TqplE4lAS_8ibglZA7iDFx7k2o/edit?usp=sharing

from gagneurlab.gagneur_utils import get_FRASER_results_spreadsheet
from gspread_dataframe import set_with_dataframe
from gspread import WorksheetNotFound

spreadsheet = get_FRASER_results_spreadsheet()

for label, results_df in results.items():
    worksheet_name = f"{label} (auto)"
    try:
        worksheet = spreadsheet.worksheet(worksheet_name)
    except WorksheetNotFound:
        worksheet = spreadsheet.add_worksheet(worksheet_name, 1, 1)
        print("Created", worksheet.title)

    set_with_dataframe(worksheet, results_df.reset_index().fillna(''), resize=True)
    print("Updated", worksheet.title)

#%%

# use OMIM API to add OMIM info


#%%

# seqr gene lists


#%%

# RNA-seq sample metadata

#%%

# add constraint info
constraint_df = pd.read_table("https://storage.googleapis.com/gnomad-public/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz", compression="gzip")

#for c in sorted(constraint_df.columns):
#    print(c)

for k in results.keys():
    result = results[k]
    result.geneID = result.geneID.apply(lambda s: s.split(".")[0])
    break

#%%

# add gene list column
