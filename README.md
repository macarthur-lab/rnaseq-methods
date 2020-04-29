### multiqc
* [batch_0](https://macarthur-lab.github.io/rnaseq-methods/pipelines/multiqc/batch_0.html)
* [batch_1_muntoni](https://macarthur-lab.github.io/rnaseq-methods/pipelines/multiqc/batch_1_muntoni.html)
* [batch_2020_04](https://macarthur-lab.github.io/rnaseq-methods/pipelines/multiqc/batch_2020_04.html)
* combined view across all batches: [all batches](https://macarthur-lab.github.io/rnaseq-methods/pipelines/multiqc/all.html)


### pipelines

Broad GP currently delivers RNA-seq data as hg19 bams.
     
---
Then the first stage of the TGG pipeline consists of the following steps run on Terra:
  1. SamToFastq 
    ([terra](https://app.terra.bio/#workspaces/macarthurlab-rnaseq-terra/macarthurlab-rnaseq-terra/workflows/broadinstitute_gtex/samtofastq_v1-0_BETA_cfg))
    ([wdl](https://portal.firecloud.org/?return=terra#methods/broadinstitute_gtex/samtofastq_v1-0_BETA/6/wdl))
  
  2. STAR alignment
    ([terra](https://app.terra.bio/#workspaces/macarthurlab-rnaseq-terra/macarthurlab-rnaseq-terra/workflows/broadinstitute_gtex/star_v1-0_BETA_cfg))
    ([wdl](https://portal.firecloud.org/?return=terra#methods/broadinstitute_gtex/star_v1-0_BETA/7/wdl))

  3a. RNAseQC
    ([terra](https://app.terra.bio/#workspaces/macarthurlab-rnaseq-terra/macarthurlab-rnaseq-terra/workflows/broadinstitute_gtex/rnaseqc2_v1-0_BETA_cfg))
    ([wdl](https://portal.firecloud.org/?return=terra#methods/broadinstitute_gtex/rnaseqc2_v1-0_BETA/2/wdl))
  
  3b. FastQC
    ([terra](https://app.terra.bio/#workspaces/macarthurlab-rnaseq-terra/macarthurlab-rnaseq-terra/workflows/sanand/FastQC))
    ([wdl](https://portal.firecloud.org/?return=terra#methods/sanand/FastQC/1/wdl))
  
---
 
Next, to copy output files from the Terra output bucket to the TGG RNA-seq bucket, run:
```
cd rnaseq_methods/pipelines
python3 ./transfer_files_to_macarthurlab_rnaseq_bucket.py -w [workspace ID] \
    macarthurlab-rnaseq-terra [RNA-seq batch name]
```
For example:
```
python3 ./transfer_files_to_macarthurlab_rnaseq_bucket.py -w a2220985-8c41-4c55-84b6-b1a219add9bf macarthurlab-rnaseq-terra batch_0
```

Then, to update the [metadata spreadsheet](https://docs.google.com/spreadsheets/d/1S3l28tZqFmzqqwqi_BCzuIkaVFmZz9eGpGtqtH5eVoo/edit#gid=421510693) 
and add the new file paths, run 

```
cd rnaseq_methods/pipelines/sample_metadata
python3 -m pip install -r requirements.txt
```
and then run through `step1_update_data_paths_worksheet.py` and `step2_update_seqr_and_other_metadata_worksheet.py` 
interactively
(TODO convert these to scripts) 
  
---

Now that all new samples are in the  [metadata spreadsheet](https://docs.google.com/spreadsheets/d/1S3l28tZqFmzqqwqi_BCzuIkaVFmZz9eGpGtqtH5eVoo/edit#gid=421510693),
run downstream analyses - using python scripts and [hail Batch](https://hail.is/docs/batch/api.html) ([zulip](https://hail.zulipchat.com/#narrow/stream/223457-Batch-support)).

- QC
    - Impute tissue 
    - Impute sex
    - Check sample ID vs. DNA
    - Impute Ancestry (?)
    
- TGG-viewer
    - add samples to config
    - TODO: reference data (GTEx, mappability)
    - TODO: gCNV tracks
- Majiq
- Fraser
- Outrider
- Aneva
- gene lists, chess genes

