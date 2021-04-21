#!/usr/bin/env python3

import argparse
import os
import pandas as pd
import re

from tgg_viewer.junctions_track_pipelines.docker.gencode_utils import parse_introns_from_gencode_gff

p = argparse.ArgumentParser(description="This script takes a tsv file (such as Fraser results) and adds a column "
    "saying whether each junction is the gff.")
p.add_argument('-g', '--gencode-gff', help="Path of gencode .gff3 file for annotating known junctions.", required=True)
p.add_argument('-c', '--chrom-column', default="seqnames")
p.add_argument('-b', '--start-column', default="start")
p.add_argument('-e', '--end-column', default="end")

g = p.add_mutually_exclusive_group(required=True)
g.add_argument("--is-0based", action="store_true")
g.add_argument("--is-1based", action="store_true")

p.add_argument('input_path', help=".tsv path")
args = p.parse_args()

input_path = args.input_path
df = pd.read_table(input_path)
print(f"Parsed {len(df)} rows and {len(df.columns)} from {input_path}")

gencode_label = re.sub(".gff3(.gz)?$", "", os.path.basename(args.gencode_gff))
gencode_introns_set = parse_introns_from_gencode_gff(args.gencode_gff)
print(f"Parsed {len(gencode_introns_set)} introns from {gencode_label}")

output_path = re.sub(".tsv(.gz)?$", "", input_path) + ".with_gencode_columns.tsv"

df.loc[:, "hidden__start"] = df[[args.start_column, args.end_column]].min(axis=1)
df.loc[:, "hidden__end"] = df[[args.start_column, args.end_column]].max(axis=1)
df.loc[:, "hidden__start_1based"] = df["hidden__start"] + 1 if args.is_0based else df["hidden__start"]
df.loc[:, "locus_1based"] = df[args.chrom_column] + ":" + df.loc[:, "hidden__start_1based"].astype("str") + "-" + df.loc[:, "hidden__end"].astype("str")
df.set_index([args.chrom_column, "hidden__start_1based", "hidden__end"], inplace=True)
df.loc[~df.index.isin(gencode_introns_set), f"isNovel (not in {gencode_label})"] = "novel"
df.fillna("", inplace=True)

print(f"Writing results to {os.path.relpath(output_path)}")
df.reset_index(inplace=True)
df.drop(columns=["hidden__start_1based", "hidden__start", "hidden__end"], inplace=True)
df.to_csv(output_path, header=True, index=False, sep="\t")

novel_junctions_counter = sum(df[f"isNovel (not in {gencode_label})"] == "novel")
print(f"{novel_junctions_counter} out of {len(df)} ({100.0*novel_junctions_counter/len(df):0.1f}%) of junctions were novel")

