import datetime
import glob
import gzip
import hail as hl
import logging
import os
import pandas as pd

from sample_metadata.rnaseq_metadata_utils import get_rnaseq_downstream_analysis_metadata_df

COUNTS_OR_TPMS = "counts"
COUNTS_OR_TPMS = "tpms"

def run(cmd):
    print(cmd)
    os.system(cmd)

#%%

hl.init(log="/dev/null")

#%%

logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s', level=logging.INFO)
logger = logging.getLogger(__name__)

if COUNTS_OR_TPMS == "tpms":
    RDG_GENE_FILE_PATH = os.path.expanduser("~/project__rnaseq/data/samples/expression/rnaseqc_tpm/*gene_tpm.gct.gz")
    GTEX_GENE_FILE_PATH = os.path.expanduser('~/project__rnaseq/data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz')
elif COUNTS_OR_TPMS == "counts":
    RDG_GENE_FILE_PATH = os.path.expanduser("~/project__rnaseq/data/samples/expression/rnaseqc_counts/*gene_reads.gct.gz")
    GTEX_GENE_FILE_PATH = os.path.expanduser('~/project__rnaseq/data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz')
else:
    raise ValueError(COUNTS_OR_TPMS)


#%%

rnaseq_sample_metadata_df = get_rnaseq_downstream_analysis_metadata_df()

#if COUNTS_OR_TPMS == "tpms":
#elif COUNTS_OR_TPMS == "counts":
#else:
#  raise ValueError(COUNTS_OR_TPMS)

print(f"Downloading {RDG_GENE_FILE_PATH} to ~/project__rnaseq/data/samples/expression/rnaseqc_tpms")
os.system(f"""cd ~/project__rnaseq/data/samples/expression/rnaseqc_tpm && ./download_RDG_samples.sh""")
print("Done downloading tpms.")

all_rdg_and_gtex_tpms_df = None
for tissue_name in ["whole_blood", "lymphocytes", "fibroblasts", "muscle"]:
    logging.info("----------------")
    logging.info(f"{tissue_name}:")

    tgg_and_gtex_metadata_df = rnaseq_sample_metadata_df[(rnaseq_sample_metadata_df["tissue"] == tissue_name)]  #  & (rnaseq_sample_metadata_df["is_GTEx_sample"] == "FALSE")
    logging.info(f"Got {len(tgg_and_gtex_metadata_df)} TGG & GTEx samples for {tissue_name}")

    # create OUTRIDER data input table
    RDG_file_paths = {}
    for file_path in glob.glob(RDG_GENE_FILE_PATH):
        sample_id = os.path.basename(file_path).replace(".gene_tpm.gct.gz", "")
        if sample_id in set(tgg_and_gtex_metadata_df.sample_id):
            RDG_file_paths[sample_id] = file_path

    RDG_sample_ids = [s for s in list(tgg_and_gtex_metadata_df.sample_id) if "GTEX" not in s]
    GTEX_sample_ids = [s for s in list(tgg_and_gtex_metadata_df.sample_id) if "GTEX" in s]

    if len(RDG_file_paths) != len(RDG_sample_ids):
        raise ValueError(f"ERROR: len(file_paths) != len(samples_df): {len(RDG_file_paths)} != {len(RDG_sample_ids)}. Missing: {set(RDG_sample_ids) - set(RDG_file_paths.keys())}" )

    print(f"Found all {len(RDG_sample_ids)} RDG samples and {len(GTEX_sample_ids)} GTEX samples for {tissue_name}")

    print(f"Reading {len(GTEX_sample_ids)} GTEX samples from {GTEX_GENE_FILE_PATH}")
    with gzip.open(os.path.expanduser(GTEX_GENE_FILE_PATH)) as f:
        next(f); next(f)
        gtex_tpms_df = pd.read_table(f, usecols=['Name', 'Description'] + list(GTEX_sample_ids))
    gtex_tpms_df.rename({'Name': 'gene_id', 'Description': 'gene_name'}, axis="columns", inplace=True)
    gtex_tpms_df.loc[:, 'gene_id'] = gtex_tpms_df['gene_id'].apply(lambda x: x.split(".")[0])
    gtex_tpms_df.set_index('gene_id', inplace=True)

    rdg_tpms_df = None
    for sample_id, file_path in RDG_file_paths.items():
        print(f"Reading {file_path}")
        with gzip.open(file_path, "rt") as f:
            next(f); next(f); next(f)
            current_gene_reads_df = pd.read_table(f, names=["gene_id", "gene_name", sample_id], usecols=['gene_id', sample_id])
            current_gene_reads_df.loc[:, 'gene_id'] = current_gene_reads_df['gene_id'].apply(lambda x: x.split(".")[0])
            current_gene_reads_df = current_gene_reads_df.set_index('gene_id')
            if rdg_tpms_df is None:
                rdg_tpms_df = current_gene_reads_df
            else:
                rdg_tpms_df = pd.merge(rdg_tpms_df, current_gene_reads_df, left_index=True, right_index=True, how="outer")

    print(len(rdg_tpms_df.index), len(set(gtex_tpms_df.index)), len(set(rdg_tpms_df.index) & set(gtex_tpms_df.index)))
    rdg_and_gtex_tpms_df = pd.merge(rdg_tpms_df, gtex_tpms_df, left_index=True, right_index=True, how="outer")
    rdg_and_gtex_tpms_df = rdg_and_gtex_tpms_df.fillna(0)

    if all_rdg_and_gtex_tpms_df is None:
        all_rdg_and_gtex_tpms_df = rdg_and_gtex_tpms_df
    else:
        all_rdg_and_gtex_tpms_df = pd.merge(all_rdg_and_gtex_tpms_df, rdg_and_gtex_tpms_df, left_index=True, right_index=True, how="outer")


#%%

timestamp = datetime.datetime.now().strftime("%Y%m%d")
tsv_output_filename = f"RDG_and_GTEX_tpms_for_all_tissues__{timestamp}.tsv.gz"
print(f"Generating {os.path.abspath(tsv_output_filename)}")
all_rdg_and_gtex_tpms_df = all_rdg_and_gtex_tpms_df.fillna(0).round(2)
all_rdg_and_gtex_tpms_df = all_rdg_and_gtex_tpms_df.reset_index()
all_rdg_and_gtex_tpms_df.to_csv(tsv_output_filename, header=True, sep="\t", index=False)
tsv_output_gs_path = f"gs://macarthurlab-rnaseq/combined_tables/{tsv_output_filename}"
hl.hadoop_copy(tsv_output_filename, tsv_output_gs_path)
print(f"Wrote {len(all_rdg_and_gtex_tpms_df)} genes x {len(all_rdg_and_gtex_tpms_df.columns)} samples to {tsv_output_gs_path}")

#%%
tsv_output_filename = f"RDG_and_GTEX_tpms_for_all_tissues__with_metadata_pivoted__{timestamp}.tsv.gz"
print(f"Generating {os.path.abspath(tsv_output_filename)}")
pivoted_df = all_rdg_and_gtex_tpms_df.set_index("gene_id").stack().to_frame("tpm")
# create pivoted df with columns "gene_id", "sample_id", "tpm"
pivoted_df = pivoted_df.reset_index().rename(columns={"level_1": "sample_id"})

meta_subset_df = rnaseq_sample_metadata_df[
    ["sample_id", "project", "read_length", "sex", "stranded", "tissue", "tissue_detail"]
].set_index("sample_id")
pivoted_df = pivoted_df.set_index("sample_id").join(meta_subset_df, how="left")
pivoted_df.reset_index().to_csv(tsv_output_filename, header=True, sep="\t", index=False)
tsv_output_gs_path = f"gs://macarthurlab-rnaseq/combined_tables/{tsv_output_filename}"
hl.hadoop_copy(tsv_output_filename, tsv_output_gs_path)
print(f"Wrote {len(all_rdg_and_gtex_tpms_df)} genes x {len(all_rdg_and_gtex_tpms_df.columns)} samples to {tsv_output_gs_path}")


#%%
