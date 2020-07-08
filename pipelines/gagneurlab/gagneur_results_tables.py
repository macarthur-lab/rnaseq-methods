#%%

import glob
import os
import pandas as pd
from pprint import pprint
from collections import defaultdict
#%%

# gsutil -m cp gs://macarthurlab-rnaseq/gagneur/outrider/results/*with*.tsv.gz .

all_tables = defaultdict(list)
results_tables = "/Users/weisburd/project__rnaseq/code/rnaseq_methods/pipelines/gagneurlab/OUTRIDER_results1/*.tsv.gz"
for path in [p for p in glob.glob(results_tables) if "metadata" not in p and "only_GTEX" not in p]:
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
#pd.read_table()

#%%

def read_table(path, padj_threshold=0.05):
    df = pd.read_table(path).set_index(["sampleID", "geneID"])
    df = df[df.padjust < padj_threshold]
    return df

results = {}
for label, tables in all_tables.items():
    if len(tables) > 2:
        raise ValueError(f"More than 2 tables found for: {label}")
    if len(tables) == 1:
        results[label] = read_table(tables[0]['path'])
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

# add gene name column
constraint_df = pd.read_table("https://storage.googleapis.com/gnomad-public/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz", compression="gzip")
# load constraint file

#for c in sorted(constraint_df.columns):
#    print(c)

for k in results.keys():
    result = results[k]
    result.geneID = result.geneID.apply(lambda s: s.split(".")[0])
    break

#%%

# add gene list column
