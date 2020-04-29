#!/usr/bin/env python3

import argparse
import os
import pyBigWig as bw
import re

p = argparse.ArgumentParser(description="This script takes a bigWig file and one or more intervals via -L and outputs a bigWig file filtered to this interval list")
p.add_argument('-L', '--interval', help="Only keep records contained in the chr:start-end interval(s)", action="append")
p.add_argument('-o', '--output-path', help="Output .bigWig file path")
p.add_argument('input_path', help="Input .bigWig file path")
args = p.parse_args()

print(f"Input: {args.input_path}")

suffix = ".bigWig"
if args.interval:
    if len(args.interval) <= 3:
        suffix = "." + "__".join([i.replace(",", "").replace(":", "-") for i in args.interval]) + suffix
    else: 
        suffix = ".filtered" + suffix

output_path = args.output_path or (re.sub(".bigWig$", "", args.input_path) + suffix)
print(output_path)

def parse_interval(i):
    
    try:
        chrom, start_end = i.split(":")
        start, end = map(int, start_end.replace(",", "").split("-"))
        return chrom, start, end
    except Exception as e:
        raise ValueError(f"Couldn't parse interval: {i}. {e}")
    
filter_intervals = [parse_interval(i) for i in args.interval]

with bw.open(args.input_path) as f, bw.open(output_path, "w") as f2:
    counter = 0
    f2.addHeader(list(f.chroms().items()))
    
    for chrom in set([i[0] for i in filter_intervals]):
        chroms = []
        starts = []
        ends = []
        values = []
        filter_intervals_for_chrom = [i for i in filter_intervals if i[0] == chrom]
        for interval in f.intervals(chrom):
            start_1based, end_1based, value = interval
            for filter_interval_chrom, filter_interval_start, filter_interval_end in filter_intervals_for_chrom:
                if filter_interval_chrom == chrom and int(start_1based) >= filter_interval_start and int(end_1based) <= filter_interval_end:
                    break
            else:
                continue

            counter += 1
            chroms.append(chrom)
            starts.append(start_1based)
            ends.append(end_1based)
            values.append(value)
    f2.addEntries(chroms, starts, ends, values)
    

print(f"Wrote {counter} intervals to {output_path}")
