import argparse
import collections
import json
import logging
import os
import pprint

from sample_metadata.utils import get_joined_metadata_df

logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s')
logger = logging.getLogger()
logger.setLevel(logging.INFO)

#%%

df = get_joined_metadata_df()

df.columns

#%%
batches = collections.Counter(df.star_pipeline_batch)

print(batches)


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
            "CanditateGenes (culprit,if solved) (Beryl)": "Canditate Genes (Beryl)",
        }
        description = "<table>"
        description += "\n".join([f"<tr><td><b>{COLUMN_LABELS.get(c, c)}:</b></td><td>{row[c]}</td></tr>" for c in [
            'batch_date_from_hg19_bam_header',
            'imputed tissue',
            'read length (rnaseqc)',
            'proj (seqr)',
            'analysis status (seqr)',
            'variant tags (seqr)',
            'coded phenotype (seqr)',
            'Include in manuscript? (Beryl:Probands)',
            'Phenotype (Beryl:Probands)',
            'Clinical Diagnosis (Beryl:Supp.)',
            'Data_type (Beryl:Probands)',
            'Genetic diagnosis Status (Beryl:Probands)',
            'CanditateGenes (culprit,if solved) (Beryl:Probands)',
            'Candidate  Variants (Beryl:Probands)',
        ] if row[c]])
        description += "</table>"

        rows_by_batch[batch_name].append({'name': row.sample_id, 'data': data, "description": description})
        rows_by_batch["all"].append({'name': row.sample_id, 'data': data, "description": description})

#%%

# handle walsh batch samples - switch bucket, add several non-RNA WGS samples that are in seqr
rows_by_batch['2020_08__walsh'].extend([
    {
        'name': 'WAL_OTH2411b_1086394128_D1',
        'description': 'WGS DNA sample',
        'data': [
            {'type': 'alignment', 'url': 'gs://tgg-rnaseq-walsh/WAL_OTH2411b_1086394128_D1.cram'},
            {'type': 'vcf', 'url': 'gs://tgg-rnaseq-walsh/WAL_OTH2411b_1086394128_D1.vcf.gz'},
        ],
    }, {
        'name': 'WAL_OTH2400_OTH2405_D1',
        'description': 'WGS DNA sample',
        'data': [
            {'type': 'alignment', 'url': 'gs://tgg-rnaseq-walsh/WAL_OTH2400_OTH2405_D1.cram'},
            {'type': 'vcf', 'url': 'gs://tgg-rnaseq-walsh/WAL_OTH2400_OTH2405_D1.vcf.gz'},
        ],
    }, {
        'name': 'linkage region genes',
        'description': "Unannotated genes in the chr20:1-2,135,825 linkage region - shared by Victor",
        'data': [
            {'type': 'bed', 'url': 'gs://seqr-reference-data/GRCh38/rna-seq/chr20_linkage_region_unnanotated_genes__grch38.bed.gz'},
        ],
    }
])

for d in rows_by_batch['2020_08__walsh']:
    d['name'] = d['name'].replace("WAL_", "").replace("_D1", "")
    for data in d['data']:
        data['url'] = data['url'].replace('macarthurlab-rnaseq', 'tgg-rnaseq-walsh')



# tgg-rnaseq-walsh

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
              "colorBy": "isAnnotatedJunction",
              "colorByNumReadsThreshold": 5,
              "labelWith": "uniqueReadCount",
              "bounceHeightBasedOn": "random",
              
              "hideAnnotated": false,
              "hideUnannotated": false,
                            
              "minUniquelyMappedReads": 0,
              "minTotalReads": 1,
              "maxFractionMultiMappedReads": 1,
              "minSplicedAlignmentOverhang": 0,
              
              "minSamplesWithThisJunction": 0,
              "maxSamplesWithThisJunction": 1000000,
              "minPercentSamplesWithThisJunction": 0,
              "maxPercentSamplesWithThisJunction": 100,
              "minJunctionEndsVisible": 0,
              
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
                        "name": "GTEx Muscle",
                        "description": "All splice junctions from all 803 GTEx v3 muscle samples. The junction-spanning read counts and read coverage are summed across all samples.",
                        "data": [
                            { "type": "coverage", "url": "gs://tgg-viewer/ref/GRCh38/gtex_v8/GTEX_muscle.803_samples.bigWig" },
                            { "type": "junctions", "url": "gs://tgg-viewer/ref/GRCh38/gtex_v8/GTEX_muscle.803_samples.junctions.bed.gz" }
                        ]
                    },
                    {
                        "name": "GTEx Blood",
                        "description": "All splice junctions from all 755 GTEx v3 blood samples. The junction-spanning read counts and read coverage are summed across all samples.",
                        "data": [
                            { "type": "coverage", "url": "gs://tgg-viewer/ref/GRCh38/gtex_v8/GTEX_blood.755_samples.bigWig" },
                            { "type": "junctions", "url": "gs://tgg-viewer/ref/GRCh38/gtex_v8/GTEX_blood.755_samples.junctions.bed.gz" }
                        ]
                    },
                    {
                        "name": "GTEx Fibs",
                        "description": "All splice junctions from all 504 GTEx v3 fibroblast samples. The junction-spanning read counts and read coverage are summed across all samples.",
                        "data": [
                            { "type": "coverage", "url": "gs://tgg-viewer/ref/GRCh38/gtex_v8/GTEX_fibs.504_samples.bigWig" },
                            { "type": "junctions", "url": "gs://tgg-viewer/ref/GRCh38/gtex_v8/GTEX_fibs.504_samples.junctions.bed.gz" }
                        ]
                    },
                    {
                        "name": "GTEx Lymph",
                        "description": "All splice junctions from all 174 GTEx v3 lymphocyte samples.<br/>The junction-spanning read counts and read coverage are summed across all samples.",
                        "data": [
                            { "type": "coverage", "url": "gs://tgg-viewer/ref/GRCh38/gtex_v8/GTEX_lymphocytes.174_samples.bigWig" },
                            { "type": "junctions", "url": "gs://tgg-viewer/ref/GRCh38/gtex_v8/GTEX_lymphocytes.174_samples.junctions.bed.gz" }
                        ]
                    },
                    {
                        "name": "GTEx Brain: Cortex",
                        "description": "All splice junctions from all 255 GTEx v3 cortex samples.<br/>The junction-spanning read counts and read coverage are summed across all samples.",
                        "data": [
                            { "type": "coverage", "url": "gs://tgg-viewer/ref/GRCh38/gtex_v8/GTEX_brain_cortex.255_samples.bigWig" },
                            { "type": "junctions", "url": "gs://tgg-viewer/ref/GRCh38/gtex_v8/GTEX_brain_cortex.255_samples.junctions.bed.gz" }
                        ]
                    },
                    {
                        "name": "GTEx Brain: Front. Cortex",
                        "description": "All splice junctions from all 209 GTEx v3 frontal cortex samples.<br/>The junction-spanning read counts and read coverage are summed across all samples.",
                        "data": [
                            { "type": "coverage", "url": "gs://tgg-viewer/ref/GRCh38/gtex_v8/GTEX_frontal_cortex.209_samples.bigWig" },
                            { "type": "junctions", "url": "gs://tgg-viewer/ref/GRCh38/gtex_v8/GTEX_frontal_cortex.209_samples.junctions.bed.gz" }
                        ]
                    }
                ]
            },
            {
                "categoryName": "GTEx Tracks (Norm.)",
                "rows": [
            
                    {
                        "name": "GTEx Muscle (Norm.)",
                        "description": "Highly expressed junctions from all 803 GTEx v3 muscle samples.<br/>The junction-spanning read counts are normalized to represent the average spanning read count per-sample (see formula below). Only junctions with rounded normalized spanning read count > 0 are included in this track.<br /><br />average_unique_reads_per_muscle_sample = (total_unqiue_reads_in_all_muscle_samples / number_of_muscle_samples)<br />per_sample_normalized_read_count = raw_read_count * average_unique_reads_per_muscle_sample / total_unqiue_reads_in_this_sample<br/>normalized read count for junction = sum(per_sample_normalized_read_counts) / number_of_muscle_samples<br/>",
                        "data": [
                            { "type": "coverage", "url": "gs://tgg-viewer/ref/GRCh38/gtex_v8/GTEX_muscle.803_samples.bigWig" },
                            { "type": "junctions", "url": "gs://tgg-viewer/ref/GRCh38/gtex_v8/GTEX_muscle.803_samples.normalized.junctions.bed.gz" }
                        ]
                    },
                    {
                        "name": "GTEx Blood (Norm.)",
                        "description": "Highly expressed junctions from all 755 GTEx v3 blood samples.<br/>The junction-spanning read counts are normalized to represent the average spanning read count per-sample (see formula below). Only junctions with rounded normalized spanning read count > 0 are included in this track.<br /><br />average_unique_reads_per_blood_sample = (total_unqiue_reads_in_all_blood_samples / number_of_blood_samples)<br />per_sample_normalized_read_count = raw_read_count * average_unique_reads_per_blood_sample / total_unqiue_reads_in_this_sample<br/>normalized read count for junction = sum(per_sample_normalized_read_counts) / number_of_blood_samples<br/>",
                        "data": [
                            { "type": "coverage", "url": "gs://tgg-viewer/ref/GRCh38/gtex_v8/GTEX_blood.755_samples.bigWig" },
                            { "type": "junctions", "url": "gs://tgg-viewer/ref/GRCh38/gtex_v8/GTEX_blood.755_samples.normalized.junctions.bed.gz" }
                        ]
                    },
                    {
                        "name": "GTEx Fibs (Norm.)",
                        "description": "Highly expressed junctions from all 504 GTEx v3 fibroblast samples.<br/>The junction-spanning read counts are normalized to represent the average spanning read count per-sample (see formula below). Only junctions with rounded normalized spanning read count > 0 are included in this track.<br /><br />average_unique_reads_per_fibs_sample = (total_unqiue_reads_in_all_fibs_samples / number_of_fibs_samples)<br />per_sample_normalized_read_count = raw_read_count * average_unique_reads_per_fibs_sample / total_unqiue_reads_in_this_sample<br/>normalized read count for junction = sum(per_sample_normalized_read_counts) / number_of_fibs_samples<br/>",
                        "data": [
                            { "type": "coverage", "url": "gs://tgg-viewer/ref/GRCh38/gtex_v8/GTEX_fibs.504_samples.bigWig" },
                            { "type": "junctions", "url": "gs://tgg-viewer/ref/GRCh38/gtex_v8/GTEX_fibs.504_samples.normalized.junctions.bed.gz" }
                        ]
                    },
                    {
                        "name": "GTEx Lymph (Norm.)",
                        "description": "Highly expressed junctions from all 174 GTEx v3 lymphocyte samples.<br />The junction-spanning read counts are normalized to represent the average spanning read count per-sample (see formula below). Only junctions with rounded normalized spanning read count > 0 are included in this track.<br /><br />average_unique_reads_per_lymph_sample = (total_unqiue_reads_in_all_lymph_samples / number_of_lymph_samples)<br />per_sample_normalized_read_count = raw_read_count * average_unique_reads_per_lymph_sample / total_unqiue_reads_in_this_sample<br/>normalized read count for junction = sum(per_sample_normalized_read_counts) / number_of_lymph_samples<br/>",
                        "data": [
                            { "type": "coverage", "url": "gs://tgg-viewer/ref/GRCh38/gtex_v8/GTEX_lymphocytes.174_samples.bigWig" },
                            { "type": "junctions", "url": "gs://tgg-viewer/ref/GRCh38/gtex_v8/GTEX_lymphocytes.174_samples.normalized.junctions.bed.gz" }
                        ]
                    },
                    {
                        "name": "GTEx Brain: Cortex (Norm.)",
                        "description": "Highly expressed junctions from all 255 GTEx v3 brain cortex samples.<br />The junction-spanning read counts are normalized to represent the average spanning read count per-sample (see formula below). Only junctions with rounded normalized spanning read count > 0 are included in this track.<br /><br />average_unique_reads_per_cortex_sample = (total_unqiue_reads_in_all_cortex_samples / number_of_cortex_samples)<br />per_sample_normalized_read_count = raw_read_count * average_unique_reads_per_cortex_sample / total_unqiue_reads_in_this_sample<br/>normalized read count for junction = sum(per_sample_normalized_read_counts) / number_of_cortex_samples<br/>",
                        "data": [
                            { "type": "coverage", "url": "gs://tgg-viewer/ref/GRCh38/gtex_v8/GTEX_brain_cortex.255_samples.bigWig" },
                            { "type": "junctions", "url": "gs://tgg-viewer/ref/GRCh38/gtex_v8/GTEX_brain_cortex.255_samples.normalized.junctions.bed.gz" }
                        ]
                    },
                    {
                        "name": "GTEx Brain: Front. Cortex (Norm.)",
                        "description": "Highly expressed junctions from all 209 GTEx v3 brain frontal cortex samples.<br />The junction-spanning read counts are normalized to represent the average spanning read count per-sample (see formula below). Only junctions with rounded normalized spanning read count > 0 are included in this track.<br /><br />average_unique_reads_per_cortex_sample = (total_unqiue_reads_in_all_cortex_samples / number_of_cortex_samples)<br />per_sample_normalized_read_count = raw_read_count * average_unique_reads_per_cortex_sample / total_unqiue_reads_in_this_sample<br/>normalized read count for junction = sum(per_sample_normalized_read_counts) / number_of_cortex_samples<br/>",
                        "data": [
                            { "type": "coverage", "url": "gs://tgg-viewer/ref/GRCh38/gtex_v8/GTEX_frontal_cortex.209_samples.bigWig" },
                            { "type": "junctions", "url": "gs://tgg-viewer/ref/GRCh38/gtex_v8/GTEX_frontal_cortex.209_samples.normalized.junctions.bed.gz" }
                        ]
                    }
                ]
            },
            {
                "categoryName": "Phenotype Tracks",
                "rows": [
                    {
                        "name": "Haploinsufficiency Genes",
                        "description": "ClinGen dosage sensitivity curation tracks from https://clinicalgenome.org/working-groups/dosage-sensitivity-curation",
                        "data": [
                            { "type": "bed", "url": "gs://tgg-viewer/ref/GRCh38/clingen/ClinGen_haploinsufficiency_gene_GRCh38.sorted.bed.gz" }
                        ]
                    },
                    {
                        "name": "Triploinsufficiency Genes",
                        "description": "ClinGen dosage sensitivity curation tracks from https://clinicalgenome.org/working-groups/dosage-sensitivity-curation",
                        "data": [
                            { "type": "bed", "url": "gs://tgg-viewer/ref/GRCh38/clingen/ClinGen_triplosensitivity_gene_GRCh38.sorted.bed.gz" }
                        ]
                    },
                    {
                        "name": "Recurrent CNVs v1.1",
                        "description": "ClinGen dosage sensitivity curation tracks from https://clinicalgenome.org/working-groups/dosage-sensitivity-curation",
                        "data": [
                            { "type": "bed", "url": "gs://tgg-viewer/ref/GRCh38/clingen/ClinGen_recurrent_CNV.V1.1.sorted.bed.gz" }
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

#%%
#%%
