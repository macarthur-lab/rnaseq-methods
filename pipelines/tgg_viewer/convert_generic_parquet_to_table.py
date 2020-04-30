import argparse
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("path", help="Path of .parquet file")
    return p.parse_args()


def main():
    args = parse_args()

    df = pd.read_parquet(args.path).reset_index()
    output_prefix = args.path.replace(".parquet", "")
    df.to_csv(f"{output_prefix}.with_header.tab", header=True, sep="\t", index=False)
    df.to_csv(f"{output_prefix}.tab", header=False, sep="\t", index=False)

    print(f"Converted {args.path} out {output_prefix}.tab")


if __name__ == "__main__":
    main()
