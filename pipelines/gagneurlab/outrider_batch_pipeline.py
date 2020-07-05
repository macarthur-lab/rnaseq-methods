
import argparse
import datetime
import hail as hl  # used for hadoop file utils
import hailtop.batch as hb
import hashlib
import logging
import os
import pandas as pd
import sys

from batch import batch_utils

#%%

OUTRIDER_BATCHES = {
    "muscle_F_101bp_without_GTEX": {"tissue": "muscle", "sex": "F", "include_GTEX_samples": "FALSE", "samples": ['B15-15_1_1', 'B15-28_1_1', 'BEG_14-4_T65', 'BON_B09-27-1_1', 'BON_B12-33-2_1', 'BON_B12-76-1_2', 'BON_B13-55_1_2', 'BON_B14-163-1_2', 'BON_B14-71-2_1', 'BON_B14-75-1_1', 'BON_B15-118_1', 'BON_B16-75-1_2', 'BON_UC473_1', 'CLA_214DF_AB_2', 'CLA_329FK_RR_2', 'CLA_62R_CaM_3', 'HK069-0177_2', 'INMR_HZ_401', 'K1157-1-4', 'LIA_EDW01_1', 'MAN_1438-01-M1', 'MBEL028_001_3', 'MTEH041_001_2', 'OUN_HK116_001', 'OUN_HK137_001', 'RGP_248_3', 'RGP_54_3_2', 'RGP_800_3_2'], },
    "muscle_F_76bp_without_GTEX": {"tissue": "muscle", "sex": "F", "include_GTEX_samples": "FALSE", "samples": ['146BO_JB_M1', '149BP_AB_M1', '164BW_KC_M1', '247DT_SH_M1', '252DX_DC_M1', '26I_SK_M1', '361AL_CC_M1', '37L_NA_M1', '65T_CR_M1', '9C_DH_M1', 'B09-24-1RNA_UNKNOWN', 'B10-02-1RNA', 'B12-21-1M', 'B13-15RNA', 'B14-117-1RNA-2', 'B14-130-1M', 'B14-70-1RNA', 'B14-78-1-U', 'CLA_79Z_CP_2', 'M_0146-01-H1', 'T238', 'T850', 'UC316-1M'], },
    # removed 'CLA_62R_CaM_2',
    "muscle_M_101bp_without_GTEX": {"tissue": "muscle", "sex": "M", "include_GTEX_samples": "FALSE", "samples": ['193CP_LS', '92AI_SM', 'B16-48_1_1', 'BB0280_CH_AffF_2', 'BEG_1025-1_T999', 'BEG_1078-1_T1071', 'BEG_1230-1_T1227', 'BEG_1438-1_T1339', 'BEG_851-1_T840', 'BEG_887-1_T1240', 'BEG_916-1_T916', 'BON_B12-74-1_1', 'BON_B14-20_1', 'BON_B14-60-1_2', 'BON_B15-76-2_2', 'BON_B16-19_1', 'BON_B16-22_1', 'BON_B16-53-1_1', 'BON_B16_30_1_1', 'BON_B16_80_1_1', 'BON_B17-23_1', 'BON_B18-24_1', 'BON_B18-25_1', 'BON_B18-54_1', 'BON_UC219-1_1', 'CLA_338FT_DM_2', 'HK018_0047_2', 'HK072-001_2', 'ICCV_458_10CC06258_02', 'INMR_GB_408', 'INMR_GI_456', 'INMR_HW_389', 'INMR_IG_561', 'INMR_IK_581', 'LIA_MAS02_2', 'LIA_TIS03_2', 'MAN_1001_01_M1_D1', 'MBEL028_002_1', 'MBRU030_2', 'MCOP008_001_2', 'MESP021_001_2', 'MESP039_2', 'MGLA003_001_2', 'MMAD002_001_2', 'OUN_HK018_0047', 'OUN_HK047_0116', 'OUN_HK079_001', 'OUN_HK080_001', 'OUN_HK081_001', 'OUN_HK112_001', 'OUN_HK124_001', 'RGP_273_3', 'RGP_56_3_3', 'UWA_FAM7_PT_D16-1243_2', 'BEG_1111-1_T1105', 'BEG_1435-01_T1370', 'HK060-0154', 'RGP_552_3'], },
    # removed 'MAN_0063_01_03',
    "muscle_M_76bp_without_GTEX": {"tissue": "muscle", "sex": "M", "include_GTEX_samples": "FALSE", "samples": ['1179-1', '1211-1', '1258-1', '126BG_CB_M1', '153BR_JB_M1', '163BV_JE_M1', '167BZ_SP_M1', '204H_AM_M1', '205E_BD_M1', '210DB_BW_M1', '227DJ_JP_M1', '250DV_LR_M1', '251DW_SD_M1', '253DY_HA_M1', '254DZ_WP_M1', '255EA_GC_M1', '373HQ_BTG_M1', '41M_MW_M1', '46N_RL_M1', '49O_NM_M1', '81AB_MM_M1', 'B09-25-1RNA', 'B09-40-1M', 'B11-11-1M', 'B11-48-1M', 'B12-30RNA', 'B13-07-1M', 'B13-07-1RNA', 'B13-52-1M', 'B14-07RNA', 'B14-48-1RNA', 'CLA_143BN_BB_2', 'CLA_179CI_GG_2', 'CLA_180CJ_DP_2', 'NH11-441', 'NH12-1413', 'T1244', 'T892', 'TOP_MAAC031_F_Muscle1', 'UC223-1RNA', 'UC305-1M', 'UC368-1M', 'UC393-1M', 'UC84-1RNA'], },
    # removed 'CLA_388HV_XX_1',
    "fibroblasts_F_without_GTEX": {"tissue": "fibroblasts", "sex": "F", "include_GTEX_samples": "FALSE", "samples": ['B11-25-1M', 'BON_B12-71_2', 'BON_B15-125_1_2', 'BON_B16-50_1_2', 'HK006_0016', 'HK010_0026', 'HK024_0065', 'HK106-001', 'INMR_IV_615', 'RGP_7_1_2', 'RGP_7_3_2', 'RGP_7_5_2', 'VCGS_DLS_FAM1_1', 'VCGS_FAM147_459_2', 'VCGS_FAM149_465_2', 'VCGS_FAM29_90', 'VCGS_FAM2_4_2', 'VCGS_FAM49_153', 'VCGS_FAM61_190', 'VCGS_FAM67_209', 'VCGS_FAM84_260', 'HK-006', 'HK-010', 'HK-061', 'HK-070', 'HK-071-001', 'HK-087', 'HK-108', 'HK-115', 'HK-134'], },
    "fibroblasts_M_without_GTEX": {"tissue": "fibroblasts", "sex": "M", "include_GTEX_samples": "FALSE", "samples": ['BON_B14-51_1_2', 'BON_B15-26_1_1', 'BON_B15-98_1_2', 'BON_B16-57_1_2', 'BON_B17-28_1', 'HK011_0029', 'HK028_0073', 'HK044_0109', 'HK085-001', 'HK088-001_2', 'HK101-001', 'MESP014_2', 'NH12-843_Fibroblasts', 'RGP_7_2_2', 'RGP_7_4_2', 'VCGS_FAM11_32_2', 'VCGS_FAM12_36_2', 'VCGS_FAM148_462_2', 'VCGS_FAM150_468_2', 'VCGS_FAM1_1_2', 'VCGS_FAM27_84_2', 'VCGS_FAM31_97_2', 'VCGS_FAM39_123', 'VCGS_FAM3_7_2', 'VCGS_FAM41_129_2', 'VCGS_FAM42_132_2', 'VCGS_FAM4_13_2', 'VCGS_FAM52_162_3', 'VCGS_FAM58_181', 'VCGS_FAM73_227_2', 'HK-015', 'HK-022', 'HK-025', 'HK-032', 'HK-035', 'HK-047', 'HK-071-004', 'HK-073', 'HK-075', 'HK-081', 'HK-085', 'HK-095', 'HK-100', 'HK-104', 'HK-107', 'HK-114', 'HK-116', 'HK-119', 'HK-124', 'HK-126', 'HK-127', 'HK-131', 'HK-133', 'HK-139'], },
    # removed 'MAN_0063_01_02',
    "lymphocytes_M_without_GTEX": {"tissue": "lymphocytes", "sex": "M", "include_GTEX_samples": "FALSE", "samples": ['VCGS_FAM22_69_2', 'VCGS_FAM26_81_2'], },
    "lymphocytes_M_with_GTEX": {"tissue": "lymphocytes", "sex": "M", "include_GTEX_samples": "TRUE", "samples": ['VCGS_FAM22_69_2', 'VCGS_FAM26_81_2'], },

    "whole_blood_with_GTEX": {"tissue": "whole_blood", "sex": "both", "include_GTEX_samples": "TRUE", "samples": ['RGP_94_3_2', 'RGP_94_4_2', 'RGP_554_3_2'], },

    "whole_blood_without_GTEX_test": {"tissue": "whole_blood", "sex": "both", "include_GTEX_samples": "FALSE", "samples": ['RGP_94_3_2', 'RGP_94_4_2', 'RGP_554_3_2'], },

}


#%%

logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s', level=logging.INFO)
logger = logging.getLogger(__name__)

GCLOUD_USER_ACCOUNT = "weisburd@broadinstitute.org"
GCLOUD_CREDENTIALS_LOCATION = "gs://weisburd-misc/creds"
GCLOUD_PROJECT = "seqr-project"

DOCKER_IMAGE = "weisburd/gagneurlab@sha256:d1f929b223cec047207faef50f1c873ebc5c5575e29a4b0fe0bc51d5b905bbed"
# https://i12g-gagneurweb.in.tum.de/public/workshops/RNAseq_ASHG19/input_data/annotations/gencode.v29lift37.annotation.txdb
GENCODE_TXDB = "gs://seqr-bw/project__rnaseq/gencode.v26.annotation.txdb"
ALL_METADATA_TSV = "gs://seqr-bw/project__rnaseq/metadata_table_for_all_RDG_and_GTEX_samples.tsv"
ALL_COUNTS_TSV_GZ = "gs://seqr-bw/project__rnaseq/OUTRIDER_input_table_RDG_and_GTEX_counts_for_all_tissues.tsv.gz"
OUTPUT_BASE_DIR = "gs://seqr-bw/project__rnaseq/OUTRIDER_results/"


def main():
    p = batch_utils.init_arg_parser(default_cpu=4, gsa_key_file=os.path.expanduser("~/.config/gcloud/misc-270914-cb9992ec9b25.json"))
    p.add_argument("batch_name", choices=OUTRIDER_BATCHES.keys(), help="Name of RNA-seq batch to process")
    args = p.parse_args()

    if not args.force:
        hl.init(log="/dev/null", quiet=True)

    batch_name = args.batch_name

    # process samples
    with batch_utils.run_batch(args) as batch:

        # inputs:
        #   metadata file, read counts file generated by export_gagneur_metadata_table.py

        j = batch_utils.init_job(batch, None, DOCKER_IMAGE if not args.raw else None, args.cpu, args.memory, disk_size=10)
        batch_utils.switch_gcloud_auth_to_user_account(j, GCLOUD_CREDENTIALS_LOCATION, GCLOUD_USER_ACCOUNT, GCLOUD_PROJECT)
        # copy inputs
        j.command(f"""gsutil -m cp {GENCODE_TXDB} .""")
        j.command(f"""gsutil -m cp {ALL_METADATA_TSV} {ALL_COUNTS_TSV_GZ} .""")
        output_file = os.path.join(OUTPUT_BASE_DIR, f"{batch_name}.RDS")

        if not args.force and hl.hadoop_is_file(output_file):
            logger.info(f"Output file exists: {output_file} . Skipping {batch_name}...")
            return

        batch_dict = OUTRIDER_BATCHES[batch_name]
        c_vector_of_sample_names = 'c("' + '", "'.join(batch_dict['samples']) + '")'
        batch_tissue = batch_dict['tissue']
        batch_sex = batch_dict['sex']
        batch_include_GTEX_samples = batch_dict["include_GTEX_samples"]

        j.command(f"""time xvfb-run Rscript -e '

# outrider 
library(OUTRIDER)
library(annotables)
library(data.table)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(purrr)
library(ggrepel)
library(plotly)
library(stringr)
library(RColorBrewer)
library(ggsci)
library(ggplot2)
library(gtable)
library(grid)
library(gridExtra)

possibleConfounders = c("tissue", "sex", "stranded", "read_length", "batch")    # "RIN"

# input tables generated by ~/project__rnaseq/code/rnaseq_methods/pipelines/gagneurlab/metadata/export_gagneur_metadata_table.py
# batches generated by ~/project__rnaseq/code/rnaseq_methods/pipelines/gagneurlab/metadata/metadata_notebook.py

GTEX_sampleIds = c()
if ({batch_include_GTEX_samples}) {{
    if (("{batch_sex}" == "M") || ("{batch_sex}" == "F")) {{
        GTEX_sampleIds = sampleInfo[(sampleInfo$sex == "{batch_sex}") & (sampleInfo$tissue == "{batch_tissue}") & grepl("GTEX", sampleInfo$sample_id)]$sample_id
    }} else {{
        GTEX_sampleIds = sampleInfo[(sampleInfo$tissue == "{batch_tissue}") & grepl("GTEX", sampleInfo$sample_id)]$sample_id    
    }}
}}

sampleLabel = "{batch_name}_"
sampleSubset = {c_vector_of_sample_names}
sampleSubset = c(sampleSubset, GTEX_sampleIds)
print(sampleSubset)

sampleInfo = fread("{os.path.basename(ALL_METADATA_TSV)}")
geneReadCounts = fread("{os.path.basename(ALL_COUNTS_TSV_GZ)}", select=c("gene_id", sampleSubset))

sampleInfo = sampleInfo[sampleInfo$sample_id %in% sampleSubset]
if (nrow(sampleInfo) != length(sampleSubset)) {{
    print(paste("ERROR: length(sampleInfo) != length(sampleSubset):", length(sampleInfo), length(sampleSubset)))
    quit("yes")
}}

sampleInfo$read_length = as.character(sampleInfo$read_length)
geneReadCounts = geneReadCounts[!grep("ERCC", geneReadCounts$geneId),]

geneIds = geneReadCounts$gene_id
colsMiusGeneId = colnames(geneReadCounts)[!colnames(geneReadCounts) %in% c("gene_id")]
geneReadCounts = geneReadCounts[,..colsMiusGeneId]
rownames(geneReadCounts) = geneIds

cnts = as.matrix(geneReadCounts)
rownames(cnts) = geneIds
ncol(cnts)
nrow(cnts)
if (ncol(cnts) != length(sampleSubset)) {{
    print(paste("ERROR: ncol(cnts) != length(sampleSubset):", ncol(cnts), length(sampleSubset)))
    quit("yes")
}}

sampleInfo[,sampleID:=sample_id]
ods <- OutriderDataSet(countData=cnts, colData=sampleInfo)

txdb <- loadDb("{os.path.basename(GENCODE_TXDB)}")
ods <- filterExpression(ods, gtfFile=txdb, filterGenes=FALSE)   #, fpkmCutoff=100)

g = plotFPKM(ods) + theme_bw() + theme(legend.position="bottom")
ggsave(file=paste(sampleLabel, "_plotFPKM.png", sep=""), g, device="png", type="cairo")

#plotExpressedGenes(ods)

ods <- estimateSizeFactors(ods)
sortedSizeFactors = sort(sizeFactors(ods))
g = ggplot(data=NULL, aes(y=sortedSizeFactors, x=1:ncol(ods))) + 
  geom_point(color="blue", size=1) + 
  labs(x="Sample rank", y="Size factors", title="Size factor distribution") + 
  geom_label_repel(aes(label=ifelse(sortedSizeFactors > 1.5, names(sortedSizeFactors), "")), 
                   nudge_x = -35, box.padding = 0.35, point.padding = 0.5, segment.color = "grey50") +
  geom_label_repel(aes(label=ifelse(sortedSizeFactors < 0.5, names(sortedSizeFactors), "")), 
                   nudge_x = 35, box.padding   = 0.35, point.padding = 0.5, segment.color = "grey50") +
  theme_bw()

ggsave(file=paste(sampleLabel, "_sizeFactors.png", sep=""), g, type="cairo")

print(sort(sizeFactors(ods))[1:5])

print(paste(length(ods), "genes before filtering"))
ods <- ods[mcols(ods)$passedFilter,]
print(paste(length(ods), "genes after filtering"))
plotCountCorHeatmap(ods, colGroups=possibleConfounders, normalized=FALSE, device="pdf", type="cairo", nRowCluster=1, nColCluster=1, filename=paste(sampleLabel, "_plotCountCorHeatmap_before_correction.pdf", sep=""))
plotCountGeneSampleHeatmap(ods, colGroups=possibleConfounders, normalized=FALSE, device="pdf", type="cairo", filename=paste(sampleLabel, "_plotCountGeneSampleHeatmap_before_correction.pdf", sep=""))

if (length(sampleSubset) > 5) {{
    ods = findEncodingDim(ods, BPPARAM=MulticoreParam(4, progressbar=TRUE))
    g = plotEncDimSearch(ods)
    ggsave(file=paste(sampleLabel, "_plotEncDimSearch", ".png", sep=""), g, type="cairo")
    optimal_q = metadata(ods)$opt
}} else {{
    optimal_q = length(sampleSubset)
}}

# increase / descrease by 25%

q = optimal_q
original_ods = ods

ods = OUTRIDER(original_ods, verbose=TRUE, iterations=15, q=q, BPPARAM=MulticoreParam(4, progressbar=TRUE))
saveRDS(ods, paste(sampleLabel, "_ods.RDS", sep=""))

plotCountCorHeatmap(ods, colGroups=possibleConfounders, normalized=TRUE, type="cairo", nRowCluster=1, nColCluster=1, main=paste("Count correlation heatmap q=", q, sep=""), filename=paste(sampleLabel, "_plotCountCorHeatmap_after_correction.pdf", sep=""))

plotCountGeneSampleHeatmap(ods, colGroups=possibleConfounders, normalized=TRUE, main=paste("Count Gene vs Sample Heatmap q=", q, sep=""), device="pdf", type="cairo", filename=paste(sampleLabel, "_plotCountGeneSampleHeatmap_after_correction.pdf", sep=""))

res = results(ods, padjCutoff=1)
res = res[,c("sampleID", "geneID", "pValue", "padjust", "zScore", "rawcounts")][order(padjust),]
res[, "q"] = q
write.table(res, file=paste(sampleLabel, "_ods__", "q", q, "_results.tsv", sep=""), quote=FALSE, sep="\\t", row.names=FALSE)
'""")

        # Create .tar.gz for export
        #j.command(f"""gsutil -m cp *.png {OUTPUT_BASE_DIR}""")
        j.command(f"""gsutil -m cp  *.tsv *.pdf *.png *.RDS {OUTPUT_BASE_DIR}""")

        logger.info(f"Output: {output_file}")


if __name__ == "__main__":
    main()

