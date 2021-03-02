import datetime
import glob
import gzip
import logging
import os
import pandas as pd

from sample_metadata.rnaseq_metadata_utils import get_rnaseq_downstream_analysis_metadata_df

logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s', level=logging.INFO)
logger = logging.getLogger(__name__)

RDG_GENE_TPMS_PATH = os.path.expanduser("~/project__rnaseq/data/samples/expression/rnaseqc_tpm/*gene_tpm.gct.gz")
GTEX_GENE_TPMS_PATH = os.path.expanduser('~/project__rnaseq/data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz')


def main():
    rnaseq_sample_metadata_df = get_rnaseq_downstream_analysis_metadata_df()

    print(f"Downloading {RDG_GENE_TPMS_PATH} to ~/project__rnaseq/data/samples/expression/rnaseqc_tpms")
    os.system(f"""cd ~/project__rnaseq/data/samples/expression/rnaseqc_tpms && gsutil cp -n {RDG_GENE_TPMS_PATH} .""")
    print("Done downloading tpms.")

    all_rdg_and_gtex_tpms_df = None
    for tissue_name in ["whole_blood", "lymphocytes", "fibroblasts", "muscle"]:
        logging.info("----------------")
        logging.info(f"{tissue_name}:")

        tgg_and_gtex_metadata_df = rnaseq_sample_metadata_df[(rnaseq_sample_metadata_df["tissue"] == tissue_name)]  #  & (rnaseq_sample_metadata_df["is_GTEx_sample"] == "FALSE")
        logging.info(f"Got {len(tgg_and_gtex_metadata_df)} TGG & GTEx samples for {tissue_name}")

        # create OUTRIDER data input table
        RDG_file_paths = {}
        for file_path in glob.glob(RDG_GENE_TPMS_PATH):
            sample_id = os.path.basename(file_path).replace(".gene_tpm.gct.gz", "")
            if sample_id in set(tgg_and_gtex_metadata_df.sample_id):
                RDG_file_paths[sample_id] = file_path

        RDG_sample_ids = [s for s in list(tgg_and_gtex_metadata_df.sample_id) if "GTEX" not in s]
        GTEX_sample_ids = [s for s in list(tgg_and_gtex_metadata_df.sample_id) if "GTEX" in s]

        if len(RDG_file_paths) != len(RDG_sample_ids):
            raise ValueError(f"ERROR: len(file_paths) != len(samples_df): {len(RDG_file_paths)} != {len(RDG_sample_ids)}. Missing: {set(RDG_sample_ids) - set(RDG_file_paths.keys())}" )

        print(f"Found all {len(RDG_sample_ids)} RDG samples and {len(GTEX_sample_ids)} GTEX samples for {tissue_name}")

        print(f"Reading {len(GTEX_sample_ids)} GTEX samples from {GTEX_GENE_TPMS_PATH}")
        with gzip.open(os.path.expanduser(GTEX_GENE_TPMS_PATH)) as f:
            next(f); next(f)
            gtex_tpms_df = pd.read_table(f, usecols=['Name'] + list(GTEX_sample_ids))
        gtex_tpms_df = gtex_tpms_df.rename({'Name': 'gene_id'}, axis="columns")
        gtex_tpms_df = gtex_tpms_df.set_index('gene_id')

        rdg_tpms_df = None
        for sample_id, file_path in RDG_file_paths.items():
            print(f"Reading {file_path}")
            with gzip.open(file_path, "rt") as f:
                next(f); next(f); next(f)
                current_gene_reads_df = pd.read_table(f, names=["gene_id", "name", sample_id], usecols=['gene_id', sample_id])
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

    def run(cmd):
        print(cmd)
        os.system(cmd)

    timestamp = datetime.datetime.now().strftime("%Y%m%d")
    tsv_output_path = f"RDG_and_GTEX_tpms_for_all_tissues__{timestamp}.tsv.gz"
    print(f"Generating {tsv_output_path}")
    all_rdg_and_gtex_tpms_df = all_rdg_and_gtex_tpms_df.fillna(0).round(2)
    all_rdg_and_gtex_tpms_df = all_rdg_and_gtex_tpms_df.reset_index()
    all_rdg_and_gtex_tpms_df['gene_id'] = all_rdg_and_gtex_tpms_df['gene_id'].apply(lambda x: x.split(".")[0])
    all_rdg_and_gtex_tpms_df.to_csv(tsv_output_path, header=True, sep="\t", index=False)
    print(f"Wrote {len(all_rdg_and_gtex_tpms_df)} genes x {len(all_rdg_and_gtex_tpms_df.columns)} samples to {tsv_output_path}")
    #run(f"gsutil cp {tsv_output_path} {CLOUD_STORAGE_BASE_DIR}/")

    tsv_output_path = f"RDG_and_GTEX_tpms_for_all_tissues__with_metadata_pivoted__{timestamp}.tsv.gz"
    print(f"Generating {tsv_output_path}")
    pivoted_df = all_rdg_and_gtex_tpms_df.set_index("gene_id").stack().to_frame("tpm")
    # create pivoted df with columns "gene_id", "sample_id", "tpm"
    pivoted_df = pivoted_df.reset_index().rename(columns={"level_1": "sample_id"})

    meta_subset_df = rnaseq_sample_metadata_df[
        ["sample_id", "project", "read_length", "sex", "stranded", "tissue", "tissue_detail"]
    ].set_index("sample_id")
    pivoted_df = pivoted_df.set_index("sample_id").join(meta_subset_df, how="left")
    pivoted_df.reset_index().to_csv(tsv_output_path, header=True, sep="\t", index=False)


if __name__ == "__main__":
    main()
