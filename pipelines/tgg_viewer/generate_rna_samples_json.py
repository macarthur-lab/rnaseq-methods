import argparse
import collections
import logging
import os
import pprint

from sample_metadata.utils import get_joined_metadata_df

logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s')
logger = logging.getLogger()
logger.setLevel(logging.INFO)

#%%
def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("-b", "--batch-name", help="optional batch name")
    args = p.parse_args()

    return args


def main():
    args = parse_args()
    logger.info("Args:\n" + pprint.pformat(args.__dict__))

    df = get_joined_metadata_df()

    logger.info("Done")


if __name__ == "__main__":
    #main()
    pass

#%%


df = get_joined_metadata_df()

df.columns

#%%
batches = collections.Counter(df.star_pipeline_batch)

print(batches)

#%%

import json
import os

#%%
os.chdir("/Users/weisburd/project__rnaseq/code/rnaseq_methods/pipelines/tgg_viewer")
print(os.getcwd())

#%%

#row.star_pipeline_batch

rows = []
for _, row in df.iterrows():
    if not row.sample_id:
        continue

    data = []
    if row.junctions_bed:
        data.append({ 'type': 'junctions', 'url': row.junctions_bed })
    if row.coverage_bigwig:
        data.append({ 'type': 'coverage', 'url': row.coverage_bigwig })
    if row.grch38_vcf and row.grch38_vcf.strip():
        data.append({ 'type': 'vcf', 'url': row.grch38_vcf })
    if row.star_bam:
        data.append({ 'type': 'alignment', 'url': row.star_bam })

    if data:
        rows.append({ 'name': row.sample_id, 'data': data })

#%%
# output settings

rnaseq_sample_rows = json.dumps(rows)

settings_json = """
{
    "genome": "hg38",
    "locus": "chr15:92,835,700-93,031,800",
    "bamOptions": {
	    "showBams": false,
        "trackHeight": 200,
        "viewAsPairs": false,
        "showSoftClips": true,
        "alignmentShading": "strand"
    },
    "sjOptions": {
        "showCoverage": true,
        "showJunctions": true,
        "trackHeight": 170,
        "colorBy": "strand",
	    "colorByNumReadsThreshold": 5,
        "thicknessBasedOn": "numUniqueReads",
        "bounceHeightBasedOn": "random",
        "labelUniqueReadCount": true,
        "labelMultiMappedReadCount": false,
        "labelTotalReadCount": false,
        "labelMotif": false,
        "labelAnnotatedJunction": false,
        "labelAnnotatedJunctionValue": " [A]",
	    "showOnlyPlusStrand": false,
	    "showOnlyMinusStrand": false,
        "hideAnnotated": false,
        "hideUnannotated": false,
        "minUniquelyMappedReads": 0,
        "minTotalReads": 1,
        "maxFractionMultiMappedReads": 1,
        "minSplicedAlignmentOverhang": 0
    },
    "vcfOptions": {
	    "showVcfs": false,
        "displayMode": "EXPANDED"
    },
    "rowsInCategories": [
        {
            "categoryName": "Reference Data",
            "rows": [
                {
                    "name": "GTEx 100 Muscle",
		            "description": "100 randomly-chosen GTEx v3 muscle samples combined by summing raw coverage values and raw splice-junction-spanning read counts across all 100 samples.",
		            "data": [
                        { "type": "coverage", "url": "gs://seqr-reference-data/GRCh38/rna-seq/GTEx_ref_data/muscle_100_GTEx_samples.bigWig" },
                        { "type": "junctions", "url": "gs://seqr-reference-data/GRCh38/rna-seq/GTEx_ref_data/muscle_100_GTEx_samples.junctions.bed.gz" }
                    ]
                },
                {
                    "name": "GTEx 100 Blood",
                    "description": "100 randomly-chosen GTEx v3 blood samples combined by summing raw coverage values and raw splice-junction-spanning read counts across all 100 samples.",
		            "data": [
                        { "type": "coverage", "url": "gs://seqr-reference-data/GRCh38/rna-seq/GTEx_ref_data/blood_100_GTEx_samples.bigWig" },
                        { "type": "junctions", "url": "gs://seqr-reference-data/GRCh38/rna-seq/GTEx_ref_data/blood_100_GTEx_samples.junctions.bed.gz" }
                    ]
                },
                {
                    "name": "GTEx 100 Fibs",
                    "description": "100 randomly-chosen GTEx v3 fibroblast samples combined by summing raw coverage values and raw splice-junction-spanning read counts across all 100 samples.",
		            "data": [
                        { "type": "coverage", "url": "gs://seqr-reference-data/GRCh38/rna-seq/GTEx_ref_data/fibs_100_GTEx_samples.bigWig" },
                        { "type": "junctions", "url": "gs://seqr-reference-data/GRCh38/rna-seq/GTEx_ref_data/fibs_100_GTEx_samples.junctions.bed.gz" }
                    ]
                },
                {
                    "name": "GTEx All Muscle",
                    "description": "All 803 GTEx v3 muscle samples combined by summing raw coverage values and raw splice-junction-spanning read counts across all samples.",
		            "data": [
                        { "type": "coverage", "url": "gs://seqr-reference-data/GRCh38/rna-seq/GTEx_ref_data/GTEX_muscle.803_samples.bigWig" },
                        { "type": "junctions", "url": "gs://seqr-reference-data/GRCh38/rna-seq/GTEx_ref_data/GTEX_muscle.803_samples.junctions.bed.gz" }
                    ]
                },
                {
                    "name": "GTEx All Blood",
                    "description": "All 755 GTEx v3 blood samples combined by summing raw coverage values and raw splice-junction-spanning read counts across all samples.",
		            "data": [
                        { "type": "coverage", "url": "gs://seqr-reference-data/GRCh38/rna-seq/GTEx_ref_data/GTEX_blood.755_samples.bigWig" },
                        { "type": "junctions", "url": "gs://seqr-reference-data/GRCh38/rna-seq/GTEx_ref_data/GTEX_blood.755_samples.junctions.bed.gz" }
                    ]
                },
                {
                    "name": "GTEx All Fibs",
                    "description": "All 504 GTEx v3 fibroblast samples combined by summing raw coverage values and raw splice-junction-spanning read counts across all samples.",
		            "data": [
                        { "type": "coverage", "url": "gs://seqr-reference-data/GRCh38/rna-seq/GTEx_ref_data/GTEX_fibs.504_samples.bigWig" },
                        { "type": "junctions", "url": "gs://seqr-reference-data/GRCh38/rna-seq/GTEx_ref_data/GTEX_fibs.504_samples.junctions.bed.gz" }
                    ]
                },
                {
                    "name": "GTEx All Muscle - Norm.",
                    "description": "All 803 GTEx v3 muscle samples combined by summing raw coverage values across all samples and also summing normalized splice-junction-spanning read counts across all samples. The normalization is done by computing the normalized read count for each junction as normalized_read_count = raw_read_count * scalar. Here scalar = average_unique_reads_per_muscle_sample / (total_unqiue_reads_in_this_sample * number_of_muscle_samples), and average_unique_reads_per_muscle_sample = (total_unqiue_reads_in_all_muscle_samples / number_of_muscle_samples)",
		            "data": [
                        { "type": "coverage", "url": "gs://seqr-reference-data/GRCh38/rna-seq/GTEx_ref_data/GTEX_muscle.803_samples.bigWig" },
                        { "type": "junctions", "url": "gs://seqr-reference-data/GRCh38/rna-seq/GTEx_ref_data/GTEX_muscle.803_samples.normalized.junctions.bed.gz" }
                    ]
                },
                {
                    "name": "GTEx All Blood - Norm.",
                    "description": "All 755 GTEx v3 blood samples combined by summing raw coverage values across all samples and also summing normalized splice-junction-spanning read counts across all samples. The normalization is done by computing the normalized read count for each junction as normalized_read_count = raw_read_count * scalar. Here scalar = average_unique_reads_per_blood_sample / (total_unqiue_reads_in_this_sample * number_of_blood_samples), and average_unique_reads_per_blood_sample = (total_unqiue_reads_in_all_blood_samples / number_of_blood_samples)",
		            "data": [
                        { "type": "coverage", "url": "gs://seqr-reference-data/GRCh38/rna-seq/GTEx_ref_data/GTEX_blood.755_samples.bigWig" },
                        { "type": "junctions", "url": "gs://seqr-reference-data/GRCh38/rna-seq/GTEx_ref_data/GTEX_blood.755_samples.normalized.junctions.bed.gz" }
                    ]
                },
                {
                    "name": "GTEx All Fibs - Norm.",
                    "description": "All 504 GTEx v3 fibroblast samples combined by summing raw coverage values across all samples and also summing normalized splice-junction-spanning read counts across all samples. The normalization is done by computing the normalized read count for each junction as normalized_read_count = raw_read_count * scalar. Here scalar = average_unique_reads_per_fibs_sample / (total_unqiue_reads_in_this_sample * number_of_fibs_samples), and average_unique_reads_per_fibs_sample = (total_unqiue_reads_in_all_fibs_samples / number_of_fibs_samples)",
		            "data": [
                        { "type": "coverage", "url": "gs://seqr-reference-data/GRCh38/rna-seq/GTEx_ref_data/GTEX_fibs.504_samples.bigWig" },
                        { "type": "junctions", "url": "gs://seqr-reference-data/GRCh38/rna-seq/GTEx_ref_data/GTEX_fibs.504_samples.normalized.junctions.bed.gz" }
                    ]
                },
                {
                    "name": "Splice AI scores - SNVs",
		            "data": [
                        { "type": "coverage", "url": "gs://seqr-reference-data/GRCh38/rna-seq/spliceai/spliceai_scores.raw.snv.hg38.all.bigWig" }
                    ]
                },
                {
                    "name": "Splice AI scores - SNVs - alt allele A",
		            "data": [
                        { "type": "coverage", "url": "gs://seqr-reference-data/GRCh38/rna-seq/spliceai/spliceai_scores.raw.snv.hg38.alt-allele-A.bigWig" }
                    ]
                },
                {
                    "name": "Splice AI scores - SNVs - alt allele C",
		            "data": [
                        { "type": "coverage", "url": "gs://seqr-reference-data/GRCh38/rna-seq/spliceai/spliceai_scores.raw.snv.hg38.alt-allele-C.bigWig" }
                    ]
                },
                {
                    "name": "Splice AI scores - SNVs - alt allele G",
		            "data": [
                        { "type": "coverage", "url": "gs://seqr-reference-data/GRCh38/rna-seq/spliceai/spliceai_scores.raw.snv.hg38.alt-allele-G.bigWig" }
                    ]
                },
                {
                    "name": "Splice AI scores - SNVs - alt allele T",
		            "data": [
                        { "type": "coverage", "url": "gs://seqr-reference-data/GRCh38/rna-seq/spliceai/spliceai_scores.raw.snv.hg38.alt-allele-T.bigWig" }
                    ]
                }
            ]
        },
        {
            "categoryName": "gCNV Batches",
            "rows": %(rnaseq_sample_rows)s
        }
    ]
}
""" % locals()

#%%
with open("settings.json", "wt") as f:
    json.dump(json.loads(settings_json), f, indent=2)

