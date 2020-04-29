import collections
import gzip
import os
import sys

#%%

junctions_by_start = collections.defaultdict(list)
junctions_by_end = collections.defaultdict(list)
for line in gzip.open(os.path.expanduser("~/project__rnaseq/data/samples/splice_junctions/9C_DH_M1.SJ.out.tab.gz"), "rt"):
    chrom, start, end = line.split("\t")[:3]
    start = int(start)
    end = int(end)
    junctions_by_start[chrom].append((start, end))
    junctions_by_end[chrom].append((start, end))

#%%

# SJ.out.tab is 1-based
exons = []
for i, chrom in enumerate(junctions_by_start):
    if "_" in chrom or chrom == "chrM":
        continue
    junctions_by_end_list = junctions_by_end[chrom]
    for i, (start, end) in enumerate(junctions_by_start[chrom]):
        for j in range(1, len(junctions_by_end_list)):
            _, current_end = junctions_by_end_list[j]

            if current_end >= start:
                _, previous_end = junctions_by_end_list[j - 1]
                junctions_by_end_list = junctions_by_end_list[max(j - 1, 0):]
                exon_size = start - (previous_end - 1)
                if exon_size > 2:
                    exons.append((exon_size, chrom, previous_end, start))

                #print("EXON: ", previous_end, start, " size", start - (previous_end - 1))
                break

#%%
exons.sort()
for size, chrom, start, end in exons[:30]:
    print(f"{chrom}:{start}-{end}  ({size} bp)")



