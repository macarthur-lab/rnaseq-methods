import argparse
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("-N", "--normalize-read-counts", action="store_true", help="whether to normalize read counts")
    p.add_argument("path", help="Path of .parquet file")
    return p.parse_args()

"""
COLUMNS_TO_OUTPUT = [
    "chrom",
    "start_1based",
    "end_1based",
    "strand",
    "intron_motif",
    "known_splice_junction",
    "unique_reads",
    "multi_mapped_reads",
    "maximum_overhang",
]
"""

def main():
    args = parse_args()

    df = pd.read_parquet(args.path).reset_index()

    ## write as tsv
    output_prefix = args.path.replace(".parquet", "")
    if args.normalize_read_counts:
        output_prefix += ".normalized_counts"

    df.to_csv(f"{output_prefix}.with_header.tab", header=True, sep="\t", index=False)
    df.to_csv(f"{output_prefix}.tab", header=False, sep="\t", index=False)


if __name__ == "__main__":
    main()
