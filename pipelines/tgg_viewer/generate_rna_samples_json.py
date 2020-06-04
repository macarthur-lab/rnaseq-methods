import argparse
import collections
import logging
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

rows_by_batch = collections.defaultdict(list)
for _, row in df.iterrows():
    if not row.sample_id:
        continue

    data = []
    if row.junctions_bed:
        data.append({'type': 'junctions', 'url': row.junctions_bed})
    if row.coverage_bigwig:
        data.append({'type': 'coverage', 'url': row.coverage_bigwig})
    if row.grch38_vcf and row.grch38_vcf.strip():
        data.append({'type': 'vcf', 'url': row.grch38_vcf})
    if row.star_bam:
        data.append({'type': 'alignment', 'url': row.star_bam})

    if data:
        batch_name = row.star_pipeline_batch.replace("batch_0", "original").replace("batch_1_", "").replace("batch_", "")

        COLUMN_LABELS = {
            "batch_date_from_hg19_bam_header": "sequencing date",
        }
        description = "<table>"
        description += "\n".join([f"<tr><td>{COLUMN_LABELS.get(c, c)}</td><td>{row[c]}</td></tr>" for c in [
            'batch_date_from_hg19_bam_header',
            'imputed tissue',
            'read length (rnaseqc)',
            'proj (seqr)',
            'analysis status (seqr)',
            'variant tags (seqr)',
            'coded phenotype (seqr)',
            'Include in manuscript? (Beryl)',
            'Phenotype (Beryl)',
            'Clinical Diagnosis (Beryl)',
            'Data_type (Beryl)',
            'Status (Beryl)',
            'CanditateGenes (culprit,if solved) (Beryl)',
            'Candidate  Variants (Beryl)',
        ]])
        description += "</table>"

        rows_by_batch[batch_name].append({'name': row.sample_id, 'data': data, "description": description})
        rows_by_batch["all"].append({'name': row.sample_id, 'data': data, "description": description})


#%%
# output settings
for batch_name, rows in rows_by_batch.items():
    rnaseq_sample_rows = json.dumps(rows)

    settings_json = """
    {
        "genome": "hg38",
        "locus": "chr21:45988674-45991233",
        "selectedRowNamesByCategoryName": {
              "Samples": [],
              "GTEx Tracks": [
                   "GTEx All Muscle - Norm.",
                   "GTEx All Blood - Norm.",
                   "GTEx All Fibs - Norm."
              ]
         },
        "selectedSamplesByCategoryNameAndRowName": {},
        "dataTypesToShow": [ "junctions", "coverage", "vcf" ],
        "bamOptions": {
            "trackHeight": 200,
            "viewAsPairs": false,
            "showSoftClips": true,
            "alignmentShading": "strand"
        },
        "sjOptions": {
              "bounceHeightBasedOn": "random",
              "colorBy": "isAnnotatedJunction",
              "colorByNumReadsThreshold": 5,
              "hideAnnotated": false,
              "hideUnannotated": false,
              "labelAnnotatedJunction": false,
              "labelAnnotatedJunctionValue": " [A]",
              "labelMotif": false,
              "labelMultiMappedReadCount": false,
              "labelTotalReadCount": false,
              "labelUniqueReadCount": true,
              "maxFractionMultiMappedReads": 1,
              "minSplicedAlignmentOverhang": 0,
              "minTotalReads": 1,
              "minUniquelyMappedReads": 0,
              "showOnlyMinusStrand": false,
              "showOnlyPlusStrand": false,
              "thicknessBasedOn": "numUniqueReads",
              "trackHeight": 170
        },
        "vcfOptions": {
            "displayMode": "EXPANDED"
        },
        "rowsInCategories": [
            {
                "categoryName": "GTEx Tracks",
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
                    }
                ]
            },
            {
                "categoryName": "Mappability Tracks",
                "rows": [
                    {
                        "name": "36-mer mappability ",
                        "description": "Mappability of 36-mers allowing for 2 mismatches. Generated using the same pipeline as the UCSC hg19 mappability tracks.",
                        "data": [
                            { "type": "coverage", "url": "gs://tgg-viewer/ref/GRCh38/mappability/GRCh38_no_alt_analysis_set_GCA_000001405.15-k36_m2.bw" }
                        ]
                    },
                    {
                        "name": "50-mer mappability ",
                        "description": "Mappability of 50-mers allowing for 2 mismatches. Generated using the same pipeline as the UCSC hg19 mappability tracks.",
                        "data": [
                            { "type": "coverage", "url": "gs://tgg-viewer/ref/GRCh38/mappability/GRCh38_no_alt_analysis_set_GCA_000001405.15-k50_m2.bw" }
                        ]
                    },
                    {
                        "name": "75-mer mappability ",
                        "description": "Mappability of 75-mers allowing for 2 mismatches. Generated using the same pipeline as the UCSC hg19 mappability tracks.",
                        "data": [
                            { "type": "coverage", "url": "gs://tgg-viewer/ref/GRCh38/mappability/GRCh38_no_alt_analysis_set_GCA_000001405.15-k75_m2.bw" }
                        ]
                    },
                    {
                        "name": "100-mer mappability ",
                        "description": "Mappability of 100-mers allowing for 2 mismatches. Generated using the same pipeline as the UCSC hg19 mappability tracks.",
                        "data": [
                            { "type": "coverage", "url": "gs://tgg-viewer/ref/GRCh38/mappability/GRCh38_no_alt_analysis_set_GCA_000001405.15-k100_m2.bw" }
                        ]
                    }
                ]
            },
            {
                "categoryName": "Samples",
                "rows": %(rnaseq_sample_rows)s
            }
        ]
    }
    """ % locals()

    output_path = "configs/"+batch_name+"_rnaseq_samples.json"
    with open(output_path, "wt") as f:
        json.dump(json.loads(settings_json), f, indent=2, sort_keys=True)
    print("Wrote to " + output_path)


#%%


x = """
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
"""
