#!/usr/bin/env bash

# This script takes a STAR splice junction file (*.SJ.out.tab) and outputs it as a .bed.gz file with 6 columns and a .tbi index which can be loaded into IGV.js

INPUT_PATH=$1
echo Input: $INPUT_PATH

OUTPUT_FORMAT='
.bed format:
column 1: chromosome
column 2: first base of the intron (0-based)
column 3: first base of the intron (1-based)
column 4: name  (overloaded to store the total number of unique + multi-mapped reads spanning the junction)
column 5: score number of uniquely-maapped reads spanning the junction
column 6: strand ("+" or "-")
'

# http://labshare.cshl.edu/shares/gingeraslab/www-data/dobin/STAR/STAR.posix/doc/STARmanual.pdf
SJ_OUT_FORMAT='
column 1: chromosome
column 2: first base of the intron (1-based)
column 3: last base of the intron (1-based)
column 4: strand (0: undefined, 1: +, 2: -)
column 5: intron motif: 0: non-canonical; 1: GT/AG, 2: CT/AC, 3: GC/AG, 4: CT/GC, 5:AT/AC, 6: GT/AT
column 6: 0: unannotated, 1: annotated (only if splice junctions database is used)
column 7: number of uniquely mapping reads crossing the junction
column 8: number of multi-mapping reads crossing the junction
column 9: maximum spliced alignment overhang
'

# https://genome.ucsc.edu/FAQ/FAQformat.html
BED_FORMAT='
chr7    127474697  127475864  Pos4  0  +
chr7    127475864  127477031  Neg1  0  -
'

OUTPUT_PATH=$(basename $INPUT_PATH | sed 's/.tab.gz//').bed.gz

# filter out junctions with 0 uniquely-mapped reads ($7 == 0)
# print chrom, start-1, end, unique+multimapped, unique, strand
echo Writing: ${OUTPUT_PATH}
gunzip -c $INPUT_PATH | awk -v OFS='\t' '{ if($7 > 0) { print }}' | awk -v OFS='\t' '{ print $1, $2 - 1, $3, $7 + $8, $7, (($4 != 2) ? "+" : "-") }' | bgzip > $OUTPUT_PATH

echo Writing: ${OUTPUT_PATH}.tbi
tabix $OUTPUT_PATH