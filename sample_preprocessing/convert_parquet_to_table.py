import argparse
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("-N", "--normalize-read-counts", action="store_true", help="whether to normalize read counts")
    p.add_argument("path", help="Path of .parquet file")
    return p.parse_args()


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


def main():
    args = parse_args()

    df = pd.read_parquet(args.path)

    num_samples = len([c for c in df.columns if c.startswith("unique_reads_")])

    ## write as tsv
    output_prefix = f"combined.{num_samples}_samples"
    if args.normalize_read_counts:
        output_prefix += ".normalized_counts"

    df[COLUMNS_TO_OUTPUT].to_csv(f"{output_prefix}.with_header.SJ.out.tab", header=True, sep="\t", index=False)
    df[COLUMNS_TO_OUTPUT].to_csv(f"{output_prefix}.SJ.out.tab", header=False, sep="\t", index=False)


if __name__ == "__main__":
    main()
