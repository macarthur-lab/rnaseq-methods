import glob
import gzip
import logging
import os
import pandas as pd

from sample_metadata.rnaseq_metadata_utils import get_rnaseq_downstream_analysis_metadata_df

logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s', level=logging.INFO)
logger = logging.getLogger(__name__)

RDG_GENE_READS_PATH = os.path.expanduser("~/project__rnaseq/data/samples/expression/rnaseqc_counts/*gene_reads.gct.gz")
GTEX_GENE_READS_PATH = os.path.expanduser('~/project__rnaseq/data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz')


def main():
    rnaseq_sample_metadata_df = get_rnaseq_downstream_analysis_metadata_df()

    print(f"Downloading {RDG_GENE_READS_PATH} to ~/project__rnaseq/data/samples/expression/rnaseqc_counts")
    os.system(f"""cd ~/project__rnaseq/data/samples/expression/rnaseqc_counts && gsutil cp -n {RDG_GENE_READS_PATH} .""")
    print("Done downloading counts.")

    all_rdg_counts_df = None
    all_rdg_and_gtex_counts_df = None
    for tissue_name in ["whole_blood", "lymphocytes", "fibroblasts", "muscle"]:
        logging.info("----------------")
        logging.info(f"{tissue_name}:")

        tgg_and_gtex_metadata_df = rnaseq_sample_metadata_df[(rnaseq_sample_metadata_df["tissue"] == tissue_name)]  #  & (rnaseq_sample_metadata_df["is_GTEx_sample"] == "FALSE")
        logging.info(f"Got {len(tgg_and_gtex_metadata_df)} TGG & GTEx samples for {tissue_name}")

        # create OUTRIDER data input table
        RDG_file_paths = {}
        for file_path in glob.glob(RDG_GENE_READS_PATH):
            sample_id = os.path.basename(file_path).replace(".gene_reads.gct.gz", "")
            if sample_id in set(tgg_and_gtex_metadata_df.sample_id):
                RDG_file_paths[sample_id] = file_path

        RDG_sample_ids = [s for s in list(tgg_and_gtex_metadata_df.sample_id) if "GTEX" not in s]
        GTEX_sample_ids = [s for s in list(tgg_and_gtex_metadata_df.sample_id) if "GTEX" in s]

        if len(RDG_file_paths) != len(RDG_sample_ids):
            raise ValueError(f"ERROR: len(file_paths) != len(samples_df): {len(RDG_file_paths)} != {len(RDG_sample_ids)}. Missing: {set(RDG_sample_ids) - set(RDG_file_paths.keys())}" )

        print(f"Found all {len(RDG_sample_ids)} RDG samples and {len(GTEX_sample_ids)} GTEX samples for {tissue_name}")

        print(f"Reading {len(GTEX_sample_ids)} GTEX samples from {GTEX_GENE_READS_PATH}")
        with gzip.open(os.path.expanduser(GTEX_GENE_READS_PATH)) as f:
            next(f); next(f)
            gtex_counts_df = pd.read_table(f, usecols=['Name'] + list(GTEX_sample_ids))
        gtex_counts_df = gtex_counts_df.rename({'Name': 'gene_id'}, axis="columns")
        gtex_counts_df = gtex_counts_df.set_index('gene_id')

        rdg_counts_df = None
        for sample_id, file_path in RDG_file_paths.items():
            print(f"Reading {file_path}")
            with gzip.open(file_path, "rt") as f:
                next(f); next(f); next(f)
                current_gene_reads_df = pd.read_table(f, names=["gene_id", "name", sample_id], usecols=['gene_id', sample_id])
                current_gene_reads_df = current_gene_reads_df.set_index('gene_id')
                if rdg_counts_df is None:
                    rdg_counts_df = current_gene_reads_df
                else:
                    rdg_counts_df = pd.merge(rdg_counts_df, current_gene_reads_df, left_index=True, right_index=True, how="outer")

        print(len(rdg_counts_df.index), len(set(gtex_counts_df.index)), len(set(rdg_counts_df.index) & set(gtex_counts_df.index)))
        rdg_and_gtex_counts_df = pd.merge(rdg_counts_df, gtex_counts_df, left_index=True, right_index=True, how="outer")
        rdg_and_gtex_counts_df = rdg_and_gtex_counts_df.fillna(0)
        #tsv_output_path = f"OUTRIDER_input_table_for_{tissue_name}.tsv"
        #rdg_and_gtex_counts_df.reset_index().to_csv(tsv_output_path, sep="\t", index=False)
        #print(f"Wrote {len(rdg_and_gtex_counts_df)} genes x {len(rdg_and_gtex_counts_df.columns)} samples to {tsv_output_path}")

        #if all_rdg_counts_df is None:
        #    all_rdg_counts_df = rdg_counts_df
        #else:
        #    all_rdg_counts_df = pd.merge(all_rdg_counts_df, rdg_counts_df, left_index=True, right_index=True, how="outer")

        if all_rdg_and_gtex_counts_df is None:
            all_rdg_and_gtex_counts_df = rdg_and_gtex_counts_df
        else:
            all_rdg_and_gtex_counts_df = pd.merge(all_rdg_and_gtex_counts_df, rdg_and_gtex_counts_df, left_index=True, right_index=True, how="outer")

    def run(cmd):
        print(cmd)
        os.system(cmd)

    #all_rdg_counts_df = all_rdg_counts_df.fillna(0)
    #tsv_output_path = f"OUTRIDER_input_table_RDG_counts_for_all_tissues.tsv"
    #all_rdg_counts_df.reset_index().to_csv(tsv_output_path, sep="\t", index=False)
    #print(f"Wrote {len(all_rdg_counts_df)} genes x {len(all_rdg_counts_df.columns)} samples to {tsv_output_path}")

    all_rdg_and_gtex_counts_df = all_rdg_and_gtex_counts_df.fillna(0)
    tsv_output_path = f"OUTRIDER_input_table_RDG_and_GTEX_counts_for_all_tissues.tsv.gz"
    all_rdg_and_gtex_counts_df.reset_index().to_csv(tsv_output_path, sep="\t", index=False)
    print(f"Wrote {len(all_rdg_and_gtex_counts_df)} genes x {len(all_rdg_and_gtex_counts_df.columns)} samples to {tsv_output_path}")
    #run(f"gsutil cp {tsv_output_path} {OUTRIDER_COUNTS_TSV_GZ}")


if __name__ == "__main__":
    main()
