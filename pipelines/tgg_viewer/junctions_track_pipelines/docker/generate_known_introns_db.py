import argparse
import pickle
from gencode_utils import parse_introns_from_gencode_gff


def _parse_args():
    p = argparse.ArgumentParser(
        description="This script takes a Gencode gff3 file and generates an intermediate .pickle file that can be passed "
                    "to convert_SJ_out_tab_to_junctions_bed.py. This is needed because in 2-pass mode STAR annotates any "
                    "junction as 'known' if it was found on the 1st pass, which represents almost all junctions.")

    p.add_argument('gencode_gff', help="Path of Gencode .gff3 file")
    args = p.parse_args()

    print(f"Input: {args.gencode_gff}")
    return args


if __name__ == "__main__":
    args = _parse_args()
    introns = parse_introns_from_gencode_gff(args.gencode_gff)
    with open(args.gencode_gff.replace(".gff3", "").replace(".gz", "") + ".pickle", "wb") as f:
        pickle.dump(set(introns), f)

