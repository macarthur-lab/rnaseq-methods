### multiqc
* [batch_0](https://macarthur-lab.github.io/rnaseq-methods/pipelines/multiqc/batch_0.html)
* [batch_1_muntoni](https://macarthur-lab.github.io/rnaseq-methods/pipelines/multiqc/batch_1_muntoni.html)
* [batch_2020_04](https://macarthur-lab.github.io/rnaseq-methods/pipelines/multiqc/batch_2020_04.html)
* [batch_2020_08](https://macarthur-lab.github.io/rnaseq-methods/pipelines/multiqc/batch_2020_08.html)
* combined view across all batches: [all batches](https://macarthur-lab.github.io/rnaseq-methods/pipelines/multiqc/all.html)


### pipelines

Broad GP currently delivers RNA-seq data as hg19 bams.
   
---
When GP delivers new samples, we 
1. manually copy the new .bam and .bai files to `gs://macarthurlab-rnaseq/[batch name]/hg19_bams/`
2. run the steps in `step1_update_data_paths_worksheet.py` to updated the google docs sheet and generate sample and sample-set metadata tables which 
can be uploaded to the [macarthurlab-rnaseq-terra](https://app.terra.bio/#workspaces/macarthurlab-rnaseq-terra/macarthurlab-rnaseq-terra/workflows) workspace for running the following workflows.  
3. update this README with the new batch      
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
python3 ./transfer_files_to_macarthurlab_rnaseq_bucket.py  macarthurlab-rnaseq-terra  batch_2020_08  \
    -w 7e14e341-78a4-4f9e-9830-df68fae4bb27 (= job id from terra Job History page) -t rnaseqc
```

---


**Metadata Spreadsheet**

To update the [metadata spreadsheet](https://docs.google.com/spreadsheets/d/1S3l28tZqFmzqqwqi_BCzuIkaVFmZz9eGpGtqtH5eVoo/edit#gid=421510693) 
and add the new file paths, run 

```
cd rnaseq_methods/pipelines/sample_metadata
python3 -m pip install -r requirements.txt
```
and then run through `step1_update_data_paths_worksheet.py` and `step2_update_seqr_and_other_metadata_worksheet.py` 
interactively
(TODO convert these to scripts) 
  
---
To generate files for the TGG-viewer, run 
1. `generate_junctions_bed_batch_pipeline.py` to generate .bed splice junction files
2. `generate_bigWig_coverage_batch_pipeline.py` to generate .bigWig coverage files 

Example:
```
 python3 ./tgg_viewer/junctions_track_pipelines/generate_junctions_bed_batch_pipeline.py -b batch_2020_08
```

Then, update the metadata paths worksheet again as described above.

---
To update the multiqc dashboard, run:

```
cd rnaseq_methods/pipelines/multiqc
python3 ./download_files_and_run_multiqc.py [batch name]
python3 ./download_files_and_run_multiqc.py all
``` 

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

---

**TGG-Viewer**

- all samples: ([tgg-viewer](http://tgg-viewer.broadinstitute.org/#locus=chr21:45988674-45991233&show=~(~'junctions~'coverage~'vcf)&selectedRows=~(GTEx*20Tracks~(~'GTEx*20All*20Muscle*20-*20Norm.~'GTEx*20All*20Blood*20-*20Norm.~'GTEx*20All*20Fibs*20-*20Norm.)~Samples~(~))&selectedSamples=~()&sjOptions=~(bounceHeightBasedOn~'random~colorBy~'isAnnotatedJunction~colorByNumReadsThreshold~5~hideAnnotated~false~hideUnannotated~false~labelAnnotatedJunction~false~labelAnnotatedJunctionValue~'*20*5bA*5d~labelMotif~false~labelMultiMappedReadCount~false~labelTotalReadCount~false~labelUniqueReadCount~true~maxFractionMultiMappedReads~1~minSplicedAlignmentOverhang~0~minTotalReads~1~minUniquelyMappedReads~0~showOnlyMinusStrand~false~showOnlyPlusStrand~false~thicknessBasedOn~'numUniqueReads~trackHeight~170)&vcfOptions=~(displayMode~'EXPANDED)&bamOptions=~(alignmentShading~'strand~showSoftClips~true~trackHeight~200~viewAsPairs~false)&gcnvOptions=~(trackHeight~200~trackMin~0~trackMax~5~autoscale~false~onlyHandleClicksForHighlightedSamples~true)&settingsUrl=~'https*3a*2f*2fraw.githubusercontent.com*2fmacarthur-lab*2frnaseq-methods*2fmaster*2fpipelines*2ftgg_viewer*2fconfigs*2fall_rnaseq_samples.json))  ([config](https://github.com/macarthur-lab/rnaseq-methods/blob/master/pipelines/tgg_viewer/configs/all_rnaseq_samples.json))
- batch_1_muntoni: ([tgg-viewer](http://tgg-viewer.broadinstitute.org/#locus=chr21:45988674-45991233&show=~(~'junctions~'coverage~'vcf)&selectedRows=~(GTEx*20Tracks~(~'GTEx*20All*20Muscle*20-*20Norm.~'GTEx*20All*20Blood*20-*20Norm.~'GTEx*20All*20Fibs*20-*20Norm.)~Samples~(~))&selectedSamples=~()&sjOptions=~(bounceHeightBasedOn~'random~colorBy~'isAnnotatedJunction~colorByNumReadsThreshold~5~hideAnnotated~false~hideUnannotated~false~labelAnnotatedJunction~false~labelAnnotatedJunctionValue~'*20*5bA*5d~labelMotif~false~labelMultiMappedReadCount~false~labelTotalReadCount~false~labelUniqueReadCount~true~maxFractionMultiMappedReads~1~minSplicedAlignmentOverhang~0~minTotalReads~1~minUniquelyMappedReads~0~showOnlyMinusStrand~false~showOnlyPlusStrand~false~thicknessBasedOn~'numUniqueReads~trackHeight~170)&vcfOptions=~(displayMode~'EXPANDED)&bamOptions=~(alignmentShading~'strand~showSoftClips~true~trackHeight~200~viewAsPairs~false)&gcnvOptions=~(trackHeight~200~trackMin~0~trackMax~5~autoscale~false~onlyHandleClicksForHighlightedSamples~true)&settingsUrl=~'https*3a*2f*2fraw.githubusercontent.com*2fmacarthur-lab*2frnaseq-methods*2fmaster*2fpipelines*2ftgg_viewer*2fconfigs*2fmuntoni_rnaseq_samples.json))  ([config](https://github.com/macarthur-lab/rnaseq-methods/blob/master/pipelines/tgg_viewer/configs/muntoni_rnaseq_samples.json))




