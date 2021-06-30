# download Gencode knownGene table from 
# https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=781208525_T8nvhQvX2wKgGYJE1KRgUZ6rIida&clade=mammal&org=Human&db=hg38&hgta_group=genes&hgta_track=knownGene&hgta_table=0&hgta_regionType=genome&position=chr1%3A11%2C102%2C837-11%2C267%2C747&hgta_outputType=bed&hgta_outFileName=

gzcat gencode_v32_knownGene.txt.gz | grep -v chrom | awk -F $'\t' 'BEGIN {OFS=FS} { print $2, $4, $5, gensub(/[.][0-9]/, "", "g", $1), "0", $3, $6, $7, "0,0,255",  $8, $9, $10 }' | grep -v _fix | grep -v _random | grep -v _alt | bedtools sort -g ~/p1/ref/GRCh38/hg38.fa.fai | bgzip > gencode_v32_knownGene.bed.gz

gzcat gencode_v32_knownGene.bed.gz | awk -F $'\t' 'BEGIN {OFS=FS} { print "0", $4, $1, $6, $2, $3, $7, $8, $10, $11, $12, "0", $4 }' | bgzip  > gencode_v32_knownGene.sorted.txt.gz

tabix -s 3 -b 5 -e 6  gencode_v32_knownGene.sorted.txt.gz 

gsutil -m cp gencode_v32_knownGene.sorted.txt.gz* gs://tgg-rnaseq/reference_tracks/



