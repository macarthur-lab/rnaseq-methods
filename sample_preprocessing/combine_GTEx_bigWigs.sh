
set -x
bigWigMerge GTEX-*.bigWig temp.bedGraph

LC_COLLATE=C sort -k1,1 -k2,2n ./temp.bedGraph > temp2.bedGraph

rm temp.bedGraph

bedGraphToBigWig temp2.bedGraph ~/p1/ref/GRCh38/hg38.fa.chromSizes merged.bigWig

rm temp2.bedGraph

