set -x

#TEST="--test"
TEST=""

for SCORE_THRESHOLD in 0.2 0.5
do 
    for GL in "gain" "loss"
    do 
      python3 convert_vcf_to_bed_and_tsv.py spliceai_scores.raw.snps_and_indels.hg38.filtered.sorted.vcf.bgz --${GL} --score-threshold ${SCORE_THRESHOLD} ${TEST} && \
      bedtools sort -i spliceai_scores.raw.snps_and_indels.hg38.filtered.sorted.score_${SCORE_THRESHOLD}.splice_${GL}.bed | bgzip > spliceai_scores.raw.snps_and_indels.hg38.filtered.sorted.score_${SCORE_THRESHOLD}.splice_${GL}.bed.gz &&\
      rm spliceai_scores.raw.snps_and_indels.hg38.filtered.sorted.score_${SCORE_THRESHOLD}.splice_${GL}.bed && \
      tabix spliceai_scores.raw.snps_and_indels.hg38.filtered.sorted.score_${SCORE_THRESHOLD}.splice_${GL}.bed.gz && \
      cp spliceai_scores.raw.snps_and_indels.hg38.filtered.sorted.score_${SCORE_THRESHOLD}.splice_${GL}.bed.gz* ~/code/igv.js-bw2/test/data/junctions/ &
    done
done

wait

gsutil -m cp spliceai_scores.raw.snps_and_indels.hg38.filtered.*splice*.bed.gz* gs://tgg-viewer/ref/GRCh38/spliceai/

