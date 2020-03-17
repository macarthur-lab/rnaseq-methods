set -x 
for allele in A C G T
do
    python3 compute_splice_ai_track.py -a $allele -r ~/p1/ref/GRCh38/hg38.fa ./spliceai_scores.raw.snv.hg38.vcf.gz >& /dev/null &
done

python3 compute_splice_ai_track.py -r ~/p1/ref/GRCh38/hg38.fa ./spliceai_scores.raw.snv.hg38.vcf.gz
