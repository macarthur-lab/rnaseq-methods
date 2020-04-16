"""
This script finds all rnaseq data files in gs://macarthurlab-rnaseq/
and updates the paths in the "Data Paths (auto)" worksheet:
https://docs.google.com/spreadsheets/d/1S3l28tZqFmzqqwqi_BCzuIkaVFmZz9eGpGtqtH5eVoo/edit#gid=216532559

sample_id, star_pipeline_batch, hg19_bam, hg19_bai, etc.
"""
#%%

from __future__ import print_function
import collections
import os
import pandas as pd
import pprint
import re
import sys

from gspread_dataframe import set_with_dataframe

from firecloud import api
from google.cloud import storage

from utils import RNASEQ_METADATA_SPREADSHEET, DATA_PATHS_WORKSHEET, RNASEQ_SAMPLE_IDS_TO_EXCLUDE

def run(cmd):
    print(cmd)
    os.system(cmd)

#%%

# Get a list of all file paths in gs://macarthurlab-rnaseq

storage_client = storage.Client()
bucket = storage_client.get_bucket('macarthurlab-rnaseq')
macarthurlab_rnaseq_bucket_blobs = bucket.list_blobs()
macarthurlab_rnaseq_bucket_file_paths = [b.public_url.replace("https://storage.googleapis.com/", "gs://") for b in macarthurlab_rnaseq_bucket_blobs]

print("Found %s paths" % len(macarthurlab_rnaseq_bucket_file_paths))


#%%

# go through all hg19 bam paths in macarthurlab_rnaseq_bucket_file_paths and initialize all_samples dict. with sample_id => {sample id, batch, hg19 bam path}
all_samples = {}

batch_sample_counters = collections.defaultdict(int)
for path in macarthurlab_rnaseq_bucket_file_paths:
    hg19_bam_path_match = re.search("macarthurlab-rnaseq/([^/]+)/hg19_bams/([^/]+).bam", path)
    if hg19_bam_path_match:
        sample_id = hg19_bam_path_match.group(2).replace(".", "-")
        if "ATYPICALMDC1A" not in path.upper() and "SIBLINGMDC1A" not in path.upper():
            sample_id = re.sub("_[TR][123]$", "", sample_id)
        sample_id = re.sub("^RP-[0-9]{0,5}_", "", sample_id)
        sample_id = re.sub("_v[1-9]_RNA_OnPrem", "", sample_id)
        sample_id = sample_id.replace("-Aligned-sortedByCoord-out", "")
        if hg19_bam_path_match.group(2) != sample_id:
            run("gsutil mv -n " + path + " " + "gs://macarthurlab-rnaseq/" + hg19_bam_path_match.group(1) + "/hg19_bams/" + sample_id + ".bam")

        if sample_id in RNASEQ_SAMPLE_IDS_TO_EXCLUDE:
            continue

        batch = hg19_bam_path_match.group(1)
        if sample_id in all_samples:
            print("ERROR: " +  sample_id + " found more than once: " + all_samples[sample_id]['star_pipeline_batch'] + ", " + batch)
            continue

        all_samples[sample_id] = {'sample_id': sample_id, 'star_pipeline_batch': batch, 'hg19_bam': path}

        batch_sample_counters[batch] += 1

pprint.pprint(dict(batch_sample_counters))

#%%

# add rest of paths in macarthurlab_rnaseq_bucket_file_paths to all_samples dict.
"""
gs://macarthurlab-rnaseq/batch_0/star/1179-1.Aligned.sortedByCoord.out.bam
gs://macarthurlab-rnaseq/batch_0/star/1179-1.Aligned.sortedByCoord.out.bam.bai
gs://macarthurlab-rnaseq/batch_0/star/1179-1.Chimeric.out.junction.gz
gs://macarthurlab-rnaseq/batch_0/star/1179-1.Log.final.out
gs://macarthurlab-rnaseq/batch_0/star/1179-1.Log.out
gs://macarthurlab-rnaseq/batch_0/star/1179-1.Log.progress.out
gs://macarthurlab-rnaseq/batch_0/star/1179-1.ReadsPerGene.out.tab.gz
gs://macarthurlab-rnaseq/batch_0/star/1179-1.SJ.out.tab.gz

gs://macarthurlab-rnaseq/batch_0/rnaseqc/1179-1.exon_reads.gct.gz
gs://macarthurlab-rnaseq/batch_0/rnaseqc/1179-1.gene_reads.gct.gz
gs://macarthurlab-rnaseq/batch_0/rnaseqc/1179-1.gene_tpm.gct.gz
gs://macarthurlab-rnaseq/batch_0/rnaseqc/1179-1.metrics.tsv
"""

for path in macarthurlab_rnaseq_bucket_file_paths:
    regexps = [
        ('hg19_bai', "macarthurlab-rnaseq/[^/]+/hg19_bams/([^/]+).bai"),

        ('star_bam', "macarthurlab-rnaseq/[^/]+/star/([^/]+).Aligned.sortedByCoord.out.bam$"),
        ('star_bai', "macarthurlab-rnaseq/[^/]+/star/([^/]+).Aligned.sortedByCoord.out.bam.bai"),
        ('star_reads_per_gene_tab', "macarthurlab-rnaseq/[^/]+/star/([^/]+).ReadsPerGene.out.tab.gz"),
        ('star_SJ_out_tab', "macarthurlab-rnaseq/[^/]+/star/([^/]+).SJ.out.tab.gz"),

        ('star_chimeric_junction', "macarthurlab-rnaseq/[^/]+/star/([^/]+).Chimeric.out.junction.gz$"),
        ('star_log_final_out', "macarthurlab-rnaseq/[^/]+/star/([^/]+).Log.final.out$"),
        ('star_log_out', "macarthurlab-rnaseq/[^/]+/star/([^/]+).Log.out$"),
        ('star_log_progress_out', "macarthurlab-rnaseq/[^/]+/star/([^/]+).Log.progress.out$"),

        ('rnaseqc_exon_reads', "macarthurlab-rnaseq/[^/]+/rnaseqc/([^/]+).exon_reads.gct.gz"),
        ('rnaseqc_gene_reads', "macarthurlab-rnaseq/[^/]+/rnaseqc/([^/]+).gene_reads.gct.gz"),
        ('rnaseqc_gene_tpm', "macarthurlab-rnaseq/[^/]+/rnaseqc/([^/]+).gene_tpm.gct.gz"),
        ('rnaseqc_metrics', "macarthurlab-rnaseq/[^/]+/rnaseqc/([^/]+).metrics.tsv"),

        ('fastqc_zip', "macarthurlab-rnaseq/[^/]+/fastqc/zip/([^/]+)_fastqc.zip"),

        ('grch38_vcf', "macarthurlab-rnaseq/[^/]+/grch38_vcfs/([^/]+).SNPs.vcf.gz$"),
        ('grch38_vcf_tbi', "macarthurlab-rnaseq/[^/]+/grch38_vcfs/([^/]+).SNPs.vcf.gz.tbi"),

        ('junctions_bed', "macarthurlab-rnaseq/[^/]+/junctions_bed_for_igv_js/([^/]+).junctions.bed.gz$"),
        ('junctions_bed_tbi', "macarthurlab-rnaseq/[^/]+/junctions_bed_for_igv_js/([^/]+).junctions.bed.gz.tbi"),

        ('coverage_bigwig', "macarthurlab-rnaseq/[^/]+/bigWig/([^/]+).bigWig"),
    ]

    for label, regexp in regexps:
        match = re.search(regexp, path)
        if not match:
            continue

        sample_id = match.group(1).replace(".", "-")
        if "ATYPICALMDC1A" not in path.upper() and "SIBLINGMDC1A" not in path.upper():
            sample_id = re.sub("_[TR][123]$", "", sample_id)
        sample_id = re.sub("^RP-[0-9]{0,5}_", "", sample_id)
        sample_id = re.sub("_v[1-9]_RNA_OnPrem", "", sample_id)
        sample_id = sample_id.replace("-Aligned-sortedByCoord-out", "")
        if match.group(1) != sample_id:
            run("gsutil mv -n " + path + " " + "/".join(path.split("/")[0:-1]) + "/" + regexp.replace("$", "").replace("([^/]+)", sample_id).split("/")[-1])

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
            print("ERROR: found more than one " + label  + " path for " +  sample_id + ":\n" +  all_samples[sample_id][label] + "\n" + path)
            continue

        all_samples[sample_id][label] = path

#%%

# use all_samples dict to init pandas DataFrame
all_samples_df = pd.DataFrame(
    columns=[
        'sample_id',
        'star_pipeline_batch',
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
        'hg19_bam',
        'hg19_bai',
        'fastqc_zip',
    ],
    data=sorted(all_samples.values(), key=lambda x: (x['star_pipeline_batch'], x['sample_id'])),
)

all_samples_df


#%%

# export the df to the  "data paths (auto)" worksheet in google docs.
set_with_dataframe(DATA_PATHS_WORKSHEET, all_samples_df.fillna(' '), resize=True)

print("Updated", RNASEQ_METADATA_SPREADSHEET.title, "/", DATA_PATHS_WORKSHEET.title)
print(DATA_PATHS_WORKSHEET.url)

#%%

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

terra_participants_table.to_csv("terra_rnaseq_samples.tsv", index=False, sep="\t")

#%%

#%%

#%%

rnaseq_workspaces = [
    ("macarthurlab-rnaseq-terra", "macarthurlab-rnaseq-terra"),
    ("macarthurlab-rnaseq-terra", "macarthurlab-rnaseq-terra - CMG_Broad_Orphan_Estonia-Onuap_RNA"),
]

header = [
    'project', 'reference_sequence_name', 'collaborator_sample_id', 'research_project', 'data_type', 'release_date', 'total_reads', 'mean_read_length', 'mean_coverage',
    'cram_path', 'crai_path', 'version', 'contamination_rate', 'chimera_rate', '20x_rate', 'library-1_mean_insert_size', 'pf_mismatch_rate', 'pf_reads',
]

total_bams = total_bais = 0
for ws_namespace, ws_name in rnaseq_workspaces:
    # returns a list of dictionaries. Each dictionary is a row in one of the metadata tables.
    entities = api.get_entities_with_type(ws_namespace, ws_name).json()
    # if the workspace has both sample, sample_set, and participant tables, the list includes rows from all of them, so filter to just the sample records
    sample_entities = [e for e in entities if e['entityType'] == "sample"]
    # for each sample, get the relevant columns including cram path, and add them to a tsv file string for upload
    tsv_string = "\t".join(["entity:sample_id", 'participant'] + header) + "\n"
    ws_total_bams = ws_total_bais = 0
    for e in sample_entities:
        sample_name = e['name']
        attr = e['attributes']
        # find which column in this workspace has the cram/bam path. Then copy it to the 'cram_path'/'crai_path' columns. I previously got this list of column names by retrieving all tables from all workspaces and looking for any column name with "cram" or "bam"
        found_bam_path = False
        for bam_key, bai_key in [
            ('bam_file', 'bai_file'),
            ('cram_or_bam_path', 'crai_or_bai_path'),
        ]:
            if attr.get(bam_key): #  and os.system("gsutil -u seqr-project ls " + attr[bam_key]) == 0:
                found_bam_path = True
                attr['hg19_bam_path'] = attr.get(bam_key)
                attr['hg19_bai_path'] = attr.get(bai_key)
                if attr.get('hg19_bam_path'):
                    break

        else:
            print("Missing bam path for: "  + sample_name)
            continue

        ws_total_bams += 1
        ws_total_bais += 1
        #participant_id = attr['participant'] if type(attr['participant']) == str else attr['participant']['entityName']
        tsv_string += "\t".join([sample_name] + list(map(str, [attr.get(k, "") for k in header]))) + "\n"
        print(tsv_string)

print("Total crams:", total_bams, "  total crais:", total_bais)

