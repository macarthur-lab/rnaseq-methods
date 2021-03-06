"""
This script finds all rnaseq data files in gs://tgg-rnaseq/
and updates the paths in the "Data Paths (auto)" worksheet:
https://docs.google.com/spreadsheets/d/1S3l28tZqFmzqqwqi_BCzuIkaVFmZz9eGpGtqtH5eVoo/edit#gid=216532559

sample_id, star_pipeline_batch, hg19_bam, hg19_bai, etc.
"""
#%%

import datetime
import collections
import hail as hl
import os
import pandas as pd
import pprint
import re
import sys

from gspread_dataframe import set_with_dataframe

from firecloud import api
from google.cloud import storage

from sample_metadata.rnaseq_metadata_utils import get_data_paths_worksheet, RNASEQ_SAMPLE_IDS_TO_EXCLUDE
from metadata.gcloud_api_utils import get_date_from_bam_header


def run(cmd):
    print(cmd)
    os.system(cmd)

hl.init(log="/dev/null")

#%%

# Get a list of all file paths in gs://tgg-rnaseq

storage_client = storage.Client()
bucket = storage_client.get_bucket('tgg-rnaseq')
tgg_rnaseq_bucket_blobs = bucket.list_blobs()
tgg_rnaseq_bucket_file_paths = [b.public_url.replace("https://storage.googleapis.com/", "gs://") for b in tgg_rnaseq_bucket_blobs]

print("Found %s paths" % len(tgg_rnaseq_bucket_file_paths))


#%%

# go through all hg19 bam paths in tgg_rnaseq_bucket_file_paths and initialize all_samples dict. with sample_id => {sample id, batch, hg19 bam path}
all_samples = {}

batch_sample_counters = collections.defaultdict(int)
for path_i, path in enumerate(tgg_rnaseq_bucket_file_paths):
    hg19_bam_path_match = re.search("tgg-rnaseq/([^/]+)/hg19_bams/([^/]+).bam$", path)
    if hg19_bam_path_match:
        sample_id = hg19_bam_path_match.group(2).replace(".", "-")
        if "ATYPICALMDC1A" not in path.upper() and "SIBLINGMDC1A" not in path.upper():
            sample_id = re.sub("_[TR][123]$", "", sample_id)
        sample_id = re.sub("^RP-[0-9]{0,5}_", "", sample_id)
        sample_id = re.sub("_v[1-9]_RNA_OnPrem", "", sample_id)
        sample_id = sample_id.replace("-Aligned-sortedByCoord-out", "")
        #print(sample_id)
        if hg19_bam_path_match.group(2) != sample_id:
            dest_path = "gs://tgg-rnaseq/" + hg19_bam_path_match.group(1) + "/hg19_bams/" + sample_id + ".bam"
            if dest_path != path and not hl.hadoop_exists(dest_path):
                #print("Would move " + path + " to " + dest_path)
                run("gsutil mv -n " + path + " " + dest_path)
            tgg_rnaseq_bucket_file_paths[path_i] = dest_path

        if sample_id in RNASEQ_SAMPLE_IDS_TO_EXCLUDE:
            continue

        batch = hg19_bam_path_match.group(1)
        if sample_id in all_samples:
            print("ERROR: " + sample_id + " found more than once in the gs://macarthurlab/*/hg19_bams/*.bam: ")
            print(all_samples[sample_id])
            print({'sample_id': sample_id, 'star_pipeline_batch': batch, 'hg19_bam': path})
            continue

        all_samples[sample_id] = {'sample_id': sample_id, 'star_pipeline_batch': batch, 'hg19_bam': path}

        batch_sample_counters[batch] += 1

pprint.pprint(dict(batch_sample_counters))


#%%

# add rest of paths in tgg_rnaseq_bucket_file_paths to all_samples dict.
"""
gs://tgg-rnaseq/batch_0/star/1179-1.Aligned.sortedByCoord.out.bam
gs://tgg-rnaseq/batch_0/star/1179-1.Aligned.sortedByCoord.out.bam.bai
gs://tgg-rnaseq/batch_0/star/1179-1.Chimeric.out.junction.gz
gs://tgg-rnaseq/batch_0/star/1179-1.Log.final.out
gs://tgg-rnaseq/batch_0/star/1179-1.Log.out
gs://tgg-rnaseq/batch_0/star/1179-1.Log.progress.out
gs://tgg-rnaseq/batch_0/star/1179-1.ReadsPerGene.out.tab.gz
gs://tgg-rnaseq/batch_0/star/1179-1.SJ.out.tab.gz

gs://tgg-rnaseq/batch_0/rnaseqc/1179-1.exon_reads.gct.gz
gs://tgg-rnaseq/batch_0/rnaseqc/1179-1.gene_reads.gct.gz
gs://tgg-rnaseq/batch_0/rnaseqc/1179-1.gene_tpm.gct.gz
gs://tgg-rnaseq/batch_0/rnaseqc/1179-1.metrics.tsv
"""

for path_i, path in enumerate(tgg_rnaseq_bucket_file_paths):
    regexps = [
        ('hg19_bai', "tgg-rnaseq/[^/]+/hg19_bams/([^/]+).bai"),

        ('star_bam', "tgg-rnaseq/[^/]+/star/([^/]+).Aligned.sortedByCoord.out.bam$"),
        ('star_bai', "tgg-rnaseq/[^/]+/star/([^/]+).Aligned.sortedByCoord.out.bam.bai"),
        ('star_reads_per_gene_tab', "tgg-rnaseq/[^/]+/star/([^/]+).ReadsPerGene.out.tab.gz"),
        ('star_SJ_out_tab', "tgg-rnaseq/[^/]+/star/([^/]+).SJ.out.tab.gz"),

        ('star_chimeric_junction', "tgg-rnaseq/[^/]+/star/([^/]+).Chimeric.out.junction.gz$"),
        ('star_log_final_out', "tgg-rnaseq/[^/]+/star/([^/]+).Log.final.out$"),
        ('star_log_out', "tgg-rnaseq/[^/]+/star/([^/]+).Log.out$"),
        ('star_log_progress_out', "tgg-rnaseq/[^/]+/star/([^/]+).Log.progress.out$"),

        ('rnaseqc_exon_reads', "tgg-rnaseq/[^/]+/rnaseqc/([^/]+).exon_reads.gct.gz"),
        ('rnaseqc_gene_reads', "tgg-rnaseq/[^/]+/rnaseqc/([^/]+).gene_reads.gct.gz"),
        ('rnaseqc_gene_tpm', "tgg-rnaseq/[^/]+/rnaseqc/([^/]+).gene_tpm.gct.gz"),
        ('rnaseqc_metrics', "tgg-rnaseq/[^/]+/rnaseqc/([^/]+).metrics.tsv"),

        ('fastqc_zip', "tgg-rnaseq/[^/]+/fastqc/zip/([^/]+)_fastqc.zip"),

        ('grch38_vcf', "tgg-rnaseq/[^/]+/grch38_vcfs/([^/]+).vcf.bgz$"),
        ('grch38_vcf_tbi', "tgg-rnaseq/[^/]+/grch38_vcfs/([^/]+).vcf.bgz.tbi"),

        ('portcullis_all', "tgg-rnaseq/[^/]+/portcullis/([^/]+).portcullis_all.junctions.tab.gz"),
        ('portcullis_filtered', "tgg-rnaseq/[^/]+/portcullis/([^/]+).portcullis_filtered.pass.junctions.tab.gz"),

        ('junctions_bed', "tgg-rnaseq/[^/]+/junctions_bed_for_igv_js/([^/]+).junctions.bed.gz$"),
        ('junctions_bed_tbi', "tgg-rnaseq/[^/]+/junctions_bed_for_igv_js/([^/]+).junctions.bed.gz.tbi"),

        ('coverage_bigwig', "tgg-rnaseq/[^/]+/bigWig/([^/]+).bigWig"),
    ]

    if "batch_all_samples" in path:
        continue

    for label, regexp in regexps:
        match = re.search(regexp, path)
        if not match:
            continue

        sample_id = match.group(1).replace(".", "-")
        if sample_id.startswith("combined-") and "junctions_bed_for_igv_js" in path:
            print(f"Skipping {path}")
            continue

        if "ATYPICALMDC1A" not in path.upper() and "SIBLINGMDC1A" not in path.upper():
            sample_id = re.sub("_[TR][123]$", "", sample_id)
        sample_id = re.sub("^RP-[0-9]{0,5}_", "", sample_id)
        sample_id = re.sub("_v[1-9]_RNA_OnPrem", "", sample_id)
        sample_id = re.sub(".SNPs", "", sample_id)
        sample_id = sample_id.replace("-Aligned-sortedByCoord-out", "")
        if match.group(1) != sample_id:
            dest_path = "/".join(path.split("/")[0:-1]) + "/" + regexp.replace("$", "").replace("([^/]+)", sample_id).split("/")[-1]
            if dest_path != path and "combined." not in path:
                #print("Would run gsutil mv -n " + path + " " + dest_path)
                run("gsutil mv -n " + path + " " + dest_path)

            tgg_rnaseq_bucket_file_paths[path_i] = dest_path


        if sample_id in RNASEQ_SAMPLE_IDS_TO_EXCLUDE:
            continue

        if sample_id not in all_samples:
            for hg19_sample_id in all_samples:
                if sample_id in hg19_sample_id:
                    sample_id = hg19_sample_id
                    #print("Found match: " + hg19_sample_id + " => " + sample_id + "  " +path)
                    break
            else:
                print("ERROR: sample_id " +  sample_id + " from " + label + " not found in all_samples dict")

            continue

        if label in all_samples[sample_id]:
            print("ERROR: found more than one " + label + " path for " + sample_id + ":\n" + all_samples[sample_id][label] + "\n" + path)
            continue

        all_samples[sample_id][label] = path

#%%

# get bam dates from header
for sample_id, record in all_samples.items():
    if record.get("hg19_bam"):
        record["batch_date_from_hg19_bam_header"] = get_date_from_bam_header(record["hg19_bam"])
        print(record["batch_date_from_hg19_bam_header"], record["hg19_bam"])

#%%

# use all_samples dict to init pandas DataFrame
all_samples_df = pd.DataFrame(
    columns=[
        'sample_id',
        'star_pipeline_batch',
        'batch_date_from_hg19_bam_header',
        'star_bam',
        'star_bai',
        'star_SJ_out_tab',
        'star_reads_per_gene_tab',
        'grch38_vcf',
        'rnaseqc_gene_reads',
        'rnaseqc_exon_reads',
        'rnaseqc_gene_tpm',
        'rnaseqc_metrics',
        'junctions_bed',
        'coverage_bigwig',
        'portcullis_filtered',
        'portcullis_all',
        'hg19_bam',
        'hg19_bai',
        'fastqc_zip',
    ],
    data=sorted(all_samples.values(), key=lambda x: (x['star_pipeline_batch'], x['batch_date_from_hg19_bam_header'], x['sample_id'])),
)

all_samples_df

#%%

print("Number of non-null data paths: ")
get_count = lambda c: sum(all_samples_df[c].str.len() > 0)
print("\n".join([f"{get_count(c)} {c}" for c in sorted(all_samples_df.columns, key=get_count, reverse=True)]))


#%%

# export the df to the  "data paths (auto)" worksheet in google docs.
set_with_dataframe(get_data_paths_worksheet(), all_samples_df.fillna(''), resize=True)

print("Updated", get_data_paths_worksheet().title)
print(get_data_paths_worksheet().url)

#%%

timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
tsv_output_path = f"~/project__rnaseq/code/rnaseq_methods/pipelines/sample_metadata/rnaseq_data_paths__{timestamp}.tsv"
all_samples_df.to_csv(tsv_output_path, sep="\t", index=False)
print(f"Wrote {len(all_samples_df)} samples to {tsv_output_path}")


#%%

os.chdir(os.path.expanduser("~/project__rnaseq/code/rnaseq_methods/pipelines/sample_metadata/"))

terra_participants_table = all_samples_df.rename({
    'sample_id': "entity:sample_id",
    'hg19_bam': "bam_file",
    'hg19_bai': "bai_file",
    'star_bam': "star_bam_file",
    'star_bai': "star_bam_index",
    'star_SJ_out_tab': "star_junctions",
    'star_reads_per_gene_tab': "star_read_counts",
    'rnaseqc_gene_tpm': "rnaseqc2_gene_tpm",
    'rnaseqc_metrics': "rnaseqc2_metrics",
    'grch38_vcf': "vcf_path",
}, axis=1)

terra_participants_table.to_csv("samples_all_batches.tsv", index=False, sep="\t")

#%%

for batch_name in set(all_samples_df.star_pipeline_batch):
    df = all_samples_df[all_samples_df.star_pipeline_batch == batch_name]
    df = df[['star_pipeline_batch', 'sample_id']].rename({
        'star_pipeline_batch': 'membership:sample_set_id',
        'sample_id': 'sample',
    }, axis=1)
    print(len(df), "samples in ", batch_name)
    df.to_csv("sample_set_" + batch_name + ".tsv", index=False, sep="\t")

#%%
