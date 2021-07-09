from collections import defaultdict
import gspread
import os
import pandas as pd
import subprocess
from google.oauth2.service_account import Credentials


"""
value_render_option: can be "FORMATTED_VALUE", "UNFORMATTED_VALUE", or "FORMULA":

    FORMATTED_VALUE: For example, if A1 is 1.23 and A2 is =A1 and formatted as currency, then A2 would return "$1.23".
    UNFORMATTED_VALUE: For example, if A1 is 1.23 and A2 is =A1 and formatted as currency, then A2 would return the number 1.23.
    FORMULA: The reply will include the formulas. For example, if A1 is 1.23 and A2 is =A1 and formatted as currency, then A2 would return "=A1".
     
https://developers.google.com/sheets/api/reference/rest/v4/ValueRenderOption
"""

VALUE_RENDER_OPTION__FORMATTED_VALUE = "FORMATTED_VALUE"
VALUE_RENDER_OPTION__UNFORMATTED_VALUE = "UNFORMATTED_VALUE"
VALUE_RENDER_OPTION__FORMULA = "FORMULA"

_GSPREAD_CLIENT = None

def get_spreasheet(spreadsheet_name):
    global _GSPREAD_CLIENT
    if _GSPREAD_CLIENT is None:
        creds = Credentials.from_service_account_file(
            os.path.expanduser('~/.config/gcloud/seqr-project-0cb2b89f436f.json'),
            scopes=[
                'https://www.googleapis.com/auth/spreadsheets',
                'https://www.googleapis.com/auth/drive.file',
                'https://www.googleapis.com/auth/drive',
            ]
        )

        _GSPREAD_CLIENT = gspread.authorize(creds)

    spreadsheet = _GSPREAD_CLIENT.open(spreadsheet_name)

    return spreadsheet

## Spreadsheet must be Shared with 733952080251-compute@developer.gserviceaccount.com
_RNASEQ_DOWNSTREAM_ANALYSIS_METADATA_WORKSHEET = None
_RNASEQ_METADATA_SPREADSHEET = None
_RNASEQ_METADATA_WORKSHEET = None
_DATA_PATHS_WORKSHEET = None
_IMPUTED_METADATA_WORKSHEET = None
_BERYLS_WORKSHEET = None
_BERYLS_WORKSHEET_2 = None
_BERYLS_WORKSHEET_3 = None

_GTEX_METADATA_SPREADSHEET = None
_GTEX_RNASEQ_SAMPLE_METADATA_WORKSHEET = None
_GTEX_WES_SAMPLE_METADATA_WORKSHEET = None
_GTEX_WGS_SAMPLE_METADATA_WORKSHEET = None
_GTEX_INDIVIDUAL_METADATA_WORKSHEET = None


def get_gtex_v8_metadata_spreadsheet():
    global _GTEX_METADATA_SPREADSHEET
    _GTEX_METADATA_SPREADSHEET = get_spreasheet("GTEx v8 metadata")
    return _GTEX_METADATA_SPREADSHEET


def get_gtex_rnaseq_sample_metadata_worksheet():
    global _GTEX_RNASEQ_SAMPLE_METADATA_WORKSHEET
    _GTEX_RNASEQ_SAMPLE_METADATA_WORKSHEET = get_gtex_v8_metadata_spreadsheet().worksheet("RNA-seq sample metadata (auto)")
    return _GTEX_RNASEQ_SAMPLE_METADATA_WORKSHEET


def get_gtex_wes_sample_metadata_worksheet():
    global _GTEX_WES_SAMPLE_METADATA_WORKSHEET
    _GTEX_WES_SAMPLE_METADATA_WORKSHEET = get_gtex_v8_metadata_spreadsheet().worksheet("WES sample metadata (auto)")
    return _GTEX_WES_SAMPLE_METADATA_WORKSHEET


def get_gtex_wgs_sample_metadata_worksheet():
    global _GTEX_WGS_SAMPLE_METADATA_WORKSHEET
    _GTEX_WGS_SAMPLE_METADATA_WORKSHEET = get_gtex_v8_metadata_spreadsheet().worksheet("WGS sample metadata (auto)")
    return _GTEX_WGS_SAMPLE_METADATA_WORKSHEET


def get_gtex_individual_metadata_worksheet():
    global _GTEX_INDIVIDUAL_METADATA_WORKSHEET
    _GTEX_INDIVIDUAL_METADATA_WORKSHEET = get_gtex_v8_metadata_spreadsheet().worksheet("individual metadata (auto)")
    return _GTEX_INDIVIDUAL_METADATA_WORKSHEET


def get_rnaseq_metadata_spreadsheet():
    global _RNASEQ_METADATA_SPREADSHEET
    _RNASEQ_METADATA_SPREADSHEET = get_spreasheet("RNA-seq metadata")
    return _RNASEQ_METADATA_SPREADSHEET


def get_rnaseq_metadata_worksheet():
    global _RNASEQ_METADATA_WORKSHEET
    _RNASEQ_METADATA_WORKSHEET = get_rnaseq_metadata_spreadsheet().worksheet("seqr info + other metadata (auto)")
    return _RNASEQ_METADATA_WORKSHEET


def get_data_paths_worksheet():
    global _DATA_PATHS_WORKSHEET
    _DATA_PATHS_WORKSHEET = get_rnaseq_metadata_spreadsheet().worksheet("data paths (auto)")
    return _DATA_PATHS_WORKSHEET


def get_rnaseq_downstream_analysis_metadata_worksheet():
    global _RNASEQ_DOWNSTREAM_ANALYSIS_METADATA_WORKSHEET
    _RNASEQ_DOWNSTREAM_ANALYSIS_METADATA_WORKSHEET = get_rnaseq_metadata_spreadsheet().worksheet("downstream analysis metadata (auto)")
    return _RNASEQ_DOWNSTREAM_ANALYSIS_METADATA_WORKSHEET


def get_imputed_metadata_worksheet():
    global _IMPUTED_METADATA_WORKSHEET
    _IMPUTED_METADATA_WORKSHEET = get_rnaseq_metadata_spreadsheet().worksheet("imputed (auto)")
    return _IMPUTED_METADATA_WORKSHEET


def get_beryls_supplementary_table_worksheet():
    global _BERYLS_WORKSHEET
    _BERYLS_WORKSHEET = get_rnaseq_metadata_spreadsheet().worksheet("Beryl's Supplementary Table 1")
    return _BERYLS_WORKSHEET


def get_beryls_rnaseq_probands_worksheet():
    global _BERYLS_WORKSHEET_2
    _BERYLS_WORKSHEET_2 = get_rnaseq_metadata_spreadsheet().worksheet("Copy of Beryl's RNAseq Probands")
    return _BERYLS_WORKSHEET_2


def get_beryls_seqr_data_worksheet():
    global _BERYLS_WORKSHEET_3
    _BERYLS_WORKSHEET_3 = get_rnaseq_metadata_spreadsheet().worksheet("Copy of Beryl's Seqr-data")
    return _BERYLS_WORKSHEET_3


def get_rnaseq_metadata_df(value_render_option=VALUE_RENDER_OPTION__FORMATTED_VALUE):
    rows = get_rnaseq_metadata_worksheet().get(value_render_option=value_render_option)
    return pd.DataFrame(data=rows[1:], columns=rows[0])

def get_rnaseq_downstream_analysis_metadata_df(value_render_option=VALUE_RENDER_OPTION__FORMATTED_VALUE):
    rows = get_rnaseq_downstream_analysis_metadata_worksheet().get(value_render_option=value_render_option)
    return pd.DataFrame(data=rows[1:], columns=rows[0])

def get_rnaseq_data_paths_df(value_render_option=VALUE_RENDER_OPTION__FORMATTED_VALUE):
    rows = get_data_paths_worksheet().get(value_render_option=value_render_option)
    return pd.DataFrame(data=rows[1:], columns=rows[0])

def get_imputed_metadata_df(value_render_option=VALUE_RENDER_OPTION__FORMATTED_VALUE):
    rows = get_imputed_metadata_worksheet().get(value_render_option=value_render_option)
    return pd.DataFrame(data=rows[1:], columns=rows[0])


def get_beryls_supplementary_table_df(value_render_option=VALUE_RENDER_OPTION__FORMATTED_VALUE):
    rows = get_beryls_supplementary_table_worksheet().get(value_render_option=value_render_option)
    return pd.DataFrame(data=rows[1:], columns=rows[0])


def get_beryls_rnaseq_probands_df(value_render_option=VALUE_RENDER_OPTION__FORMATTED_VALUE):
    rows = get_beryls_rnaseq_probands_worksheet().get(value_render_option=value_render_option)
    return pd.DataFrame(data=rows[1:], columns=rows[0])


def get_beryls_seqr_data_df(value_render_option=VALUE_RENDER_OPTION__FORMATTED_VALUE):
    rows = get_beryls_seqr_data_worksheet().get(value_render_option=value_render_option)
    return pd.DataFrame(data=rows[1:], columns=rows[0])


def get_rnaseq_metadata_joined_with_paths_df(value_render_option=VALUE_RENDER_OPTION__FORMATTED_VALUE):
    df1 = get_rnaseq_data_paths_df(value_render_option=value_render_option)
    df2 = get_rnaseq_metadata_df(value_render_option=value_render_option)
    df2 = df2[[c for c in df2.columns if c not in ("star_pipeline_batch", "batch_date_from_hg19_bam_header")]]  # remove columns that exist in both tables
    return df1.merge(df2, on="sample_id", how="left").set_index("sample_id", drop=False)


def get_rnaseqc_metrics(rnaseqc_metrics_file_path):
    output = subprocess.check_output("gsutil cat %s" % rnaseqc_metrics_file_path, shell=True, encoding="UTF-8")
    metrics_dict = {}
    for i, line in enumerate(output.rstrip().split("\n")):
        key, value = line.split("\t")
        metrics_dict[key] = value

    return metrics_dict


RNASEQ_SAMPLE_IDS_TO_EXCLUDE = {
    "VIL_17_097",
    "VIL_17_098",
    "VIL_17_099",
    "VIL_17_105",
    "VIL_17_106",
    "VIL_17_107",
    "VIL_17_110",
    "VIL_17_111",
    "VIL_17_120",
    "VIL_17_121",
    "VIL_17_141",
    "VIL_17_150",
    "VIL_17_151",
    "VIL_17_152",
    "VIL_17_163",
    "VIL_17_164",
    "VIL_18_045",
    "VIL_18_046",
    "VIL_18_061",
    "VIL_18_120",

    "HF_1",
    "HF_10",
    "HF_11",
    "HF_12",
    "HF_13",
    "HF_14",
    "HF_2",
    "HF_3",
    "HF_4",
    "HF_5",
    "HF_6",
    "HF_7",
    "HF_8",
    "HF_9",
}


def get_gtex_rnaseq_sample_metadata_df(value_render_option=VALUE_RENDER_OPTION__FORMATTED_VALUE):
    rows = get_gtex_rnaseq_sample_metadata_worksheet().get(value_render_option=value_render_option)
    return pd.DataFrame(data=rows[1:], columns=rows[0]).set_index("SAMPID", drop=False)


def get_gtex_wes_sample_metadata_df(value_render_option=VALUE_RENDER_OPTION__FORMATTED_VALUE):
    rows = get_gtex_wes_sample_metadata_worksheet().get(value_render_option=value_render_option)
    return pd.DataFrame(data=rows[1:], columns=rows[0]).set_index("SAMPID", drop=False)


def get_gtex_wgs_sample_metadata_df(value_render_option=VALUE_RENDER_OPTION__FORMATTED_VALUE):
    rows = get_gtex_wgs_sample_metadata_worksheet().get(value_render_option=value_render_option)
    return pd.DataFrame(data=rows[1:], columns=rows[0]).set_index("SAMPID", drop=False)


def get_analysis_batches():
    df = get_rnaseq_downstream_analysis_metadata_df(value_render_option=VALUE_RENDER_OPTION__FORMATTED_VALUE)
    analysis_batch_to_tissue = defaultdict(set)
    analysis_batch_to_sex = defaultdict(set)

    analysis_batches = {}
    for _, r in df.iterrows():
        analysis_batch = r["tissue"]
        if not analysis_batch:
            continue
        analysis_batch = analysis_batch.strip()
        if not analysis_batch or analysis_batch == "x":
            continue

        analysis_batch_to_tissue[analysis_batch].add(r["tissue"])
        analysis_batch_to_sex[analysis_batch].add(r["sex"])

    for analysis_batch, tissue in analysis_batch_to_tissue.items():
        if len(tissue) != 1:
            raise ValueError(f"Expected 1 tissue for {analysis_batch}. Found: {tissue}")
        tissue = next(iter(tissue))
        sex = analysis_batch_to_tissue[analysis_batch]
        if len(sex) > 1:
            sex = "both"
        else:
            sex = next(iter(sex))

        analysis_batches[analysis_batch] = {
            "tissue": tissue,
            "sex": sex,
            "samples": list(df[df["tissue"] == analysis_batch].sample_id)
        }

    # TODO fix empty values in spreadsheet "analysis batch" column
    return analysis_batches

#a = get_anaysis_batches()
#for batch in a: print(batch, set(ANALYSIS_BATCHES[batch]["samples"]) - set(a[batch]["samples"]))
"""
ANALYSIS_BATCHES = {
    "muscle_F_101bp": {"tissue": "muscle", "sex": "F", "samples": ['B15-15_1_1', 'B15-28_1_1', 'BEG_14-4_T65', 'BON_B09-27-1_1', 'BON_B12-33-2_1', 'BON_B12-76-1_2', 'BON_B13-55_1_2', 'BON_B14-163-1_2', 'BON_B14-71-2_1', 'BON_B14-75-1_1', 'BON_B15-118_1', 'BON_B16-75-1_2', 'BON_UC473_1', 'CLA_214DF_AB_2', 'CLA_329FK_RR_2', 'CLA_62R_CaM_3', 'HK069-0177_2', 'INMR_HZ_401', 'K1157-1-4', 'LIA_EDW01_1', 'MAN_1438-01-M1', 'MBEL028_001_3', 'MTEH041_001_2', 'OUN_HK116_001', 'OUN_HK137_001', 'RGP_248_3', 'RGP_54_3_2', 'RGP_800_3_2', 'BON_B09-8_1', 'BON_B19-57_1', 'RGP_800_3_3'], },
    "muscle_F_76bp": {"tissue": "muscle", "sex": "F", "samples": ['146BO_JB_M1', '149BP_AB_M1', '164BW_KC_M1', '247DT_SH_M1', '252DX_DC_M1', '26I_SK_M1', '361AL_CC_M1', '37L_NA_M1', '65T_CR_M1', '9C_DH_M1', 'B09-24-1RNA_UNKNOWN', 'B10-02-1RNA', 'B12-21-1M', 'B13-15RNA', 'B14-117-1RNA-2', 'B14-130-1M', 'B14-70-1RNA', 'B14-78-1-U', 'CLA_79Z_CP_2', 'M_0146-01-H1', 'T238', 'T850', 'UC316-1M'], },
    # removed 'CLA_62R_CaM_2',
    "muscle_M_101bp": {"tissue": "muscle", "sex": "M", "samples": ['193CP_LS', '92AI_SM', 'B16-48_1_1', 'BB0280_CH_AffF_2', 'BEG_1025-1_T999', 'BEG_1078-1_T1071', 'BEG_1230-1_T1227', 'BEG_1438-1_T1339', 'BEG_851-1_T840', 'BEG_887-1_T1240', 'BEG_916-1_T916', 'BON_B12-74-1_1', 'BON_B14-20_1', 'BON_B14-60-1_2', 'BON_B15-76-2_2', 'BON_B16-19_1', 'BON_B16-22_1', 'BON_B16-53-1_1', 'BON_B16_30_1_1', 'BON_B16_80_1_1', 'BON_B17-23_1', 'BON_B18-24_1', 'BON_B18-25_1', 'BON_B18-54_1', 'BON_UC219-1_1', 'CLA_338FT_DM_2', 'HK018_0047_2', 'HK072-001_2', 'ICCV_458_10CC06258_02', 'INMR_GB_408', 'INMR_GI_456', 'INMR_HW_389', 'INMR_IG_561', 'INMR_IK_581', 'LIA_MAS02_2', 'LIA_TIS03_2', 'MAN_1001_01_M1_D1', 'MBEL028_002_1', 'MBRU030_2', 'MCOP008_001_2', 'MESP021_001_2', 'MESP039_2', 'MGLA003_001_2', 'MMAD002_001_2', 'OUN_HK018_0047', 'OUN_HK047_0116', 'OUN_HK079_001', 'OUN_HK080_001', 'OUN_HK081_001', 'OUN_HK112_001', 'OUN_HK124_001', 'RGP_273_3', 'RGP_56_3_3', 'UWA_FAM7_PT_D16-1243_2', 'BEG_1111-1_T1105', 'BEG_1435-01_T1370', 'HK060-0154', 'RGP_552_3'], },
    # removed 'MAN_0063_01_03',
    "muscle_M_76bp": {"tissue": "muscle", "sex": "M", "samples": ['1179-1', '1211-1', '1258-1', '126BG_CB_M1', '153BR_JB_M1', '163BV_JE_M1', '167BZ_SP_M1', '204H_AM_M1', '205E_BD_M1', '210DB_BW_M1', '227DJ_JP_M1', '250DV_LR_M1', '251DW_SD_M1', '253DY_HA_M1', '254DZ_WP_M1', '255EA_GC_M1', '373HQ_BTG_M1', '41M_MW_M1', '46N_RL_M1', '49O_NM_M1', '81AB_MM_M1', 'B09-25-1RNA', 'B09-40-1M', 'B11-11-1M', 'B11-48-1M', 'B12-30RNA', 'B13-07-1M', 'B13-07-1RNA', 'B13-52-1M', 'B14-07RNA', 'B14-48-1RNA', 'CLA_143BN_BB_2', 'CLA_179CI_GG_2', 'CLA_180CJ_DP_2', 'NH11-441', 'NH12-1413', 'T1244', 'T892', 'TOP_MAAC031_F_Muscle1', 'UC223-1RNA', 'UC305-1M', 'UC368-1M', 'UC393-1M', 'UC84-1RNA'], },
    # removed 'CLA_388HV_XX_1',
    "fibroblasts_F": {"tissue": "fibroblasts", "sex": "F", "samples": ['B11-25-1M', 'BON_B12-71_2', 'BON_B15-125_1_2', 'BON_B16-50_1_2', 'HK006_0016', 'HK010_0026', 'HK024_0065', 'HK106-001', 'INMR_IV_615', 'RGP_7_1_2', 'RGP_7_3_2', 'RGP_7_5_2', 'VCGS_DLS_FAM1_1', 'VCGS_FAM147_459_2', 'VCGS_FAM149_465_2', 'VCGS_FAM29_90', 'VCGS_FAM2_4_2', 'VCGS_FAM49_153', 'VCGS_FAM61_190', 'VCGS_FAM67_209', 'VCGS_FAM84_260', 'HK-006', 'HK-010', 'HK-061', 'HK-070', 'HK-071-001', 'HK-087', 'HK-108', 'HK-115', 'HK-134', 'HK044_0110', 'BON_B19-45_CHX_1',], },
    "fibroblasts_M": {"tissue": "fibroblasts", "sex": "M", "samples": ['BON_B14-51_1_2', 'BON_B15-26_1_1', 'BON_B15-98_1_2', 'BON_B16-57_1_2', 'BON_B17-28_1', 'HK011_0029', 'HK028_0073', 'HK044_0109', 'HK085-001', 'HK088-001_2', 'HK101-001', 'MESP014_2', 'NH12-843_Fibroblasts', 'RGP_7_2_2', 'RGP_7_4_2', 'VCGS_FAM11_32_2', 'VCGS_FAM12_36_2', 'VCGS_FAM148_462_2', 'VCGS_FAM150_468_2', 'VCGS_FAM1_1_2', 'VCGS_FAM27_84_2', 'VCGS_FAM31_97_2', 'VCGS_FAM39_123', 'VCGS_FAM3_7_2', 'VCGS_FAM41_129_2', 'VCGS_FAM42_132_2', 'VCGS_FAM4_13_2', 'VCGS_FAM52_162_3', 'VCGS_FAM58_181', 'VCGS_FAM73_227_2', 'HK-015', 'HK-022', 'HK-025', 'HK-032', 'HK-035', 'HK-047', 'HK-071-004', 'HK-073', 'HK-075', 'HK-081', 'HK-085', 'HK-095', 'HK-100', 'HK-104', 'HK-107', 'HK-114', 'HK-116', 'HK-119', 'HK-124', 'HK-126', 'HK-127', 'HK-131', 'HK-133', 'HK-139', 'HK117_001', 'BON_B15-95_CHX_1', 'BON_B18-23_CHX_1', 'BON_B19-54_CHX_1', 'BON_B19-59_CHX_1',], },
    # removed 'MAN_0063_01_02',
    "lymphocytes": {"tissue": "lymphocytes", "sex": "both", "samples": ['VCGS_FAM22_69_2', 'VCGS_FAM26_81_2'], },
    "whole_blood": {"tissue": "whole_blood", "sex": "both", "samples": ['RGP_94_3_2', 'RGP_94_4_2', 'RGP_554_3_2', 'RGP_558_3', "RGP_673_3", "RGP_677_3", "RGP_677_4", ], },

    "test_batch": {"tissue": "muscle", "sex": "F", "samples": ['test1', 'test2'], },
}

SAMPLE_ID_TO_ANALYSIS_BATCH = {
    sample_id: batch_name for batch_name, d in ANALYSIS_BATCHES.items() for sample_id in d["samples"]
}

ANALYSIS_BATCHES["all_analysis_samples"] = {
    "tissue": "all_analysis_samples",
    "samples": [sample_id for batch_name in ANALYSIS_BATCHES for sample_id in ANALYSIS_BATCHES[batch_name]["samples"] if "test" not in batch_name],
    "sex": "both",
}

ANALYSIS_BATCHES["muscle"] = {
    "tissue": "muscle",
    "samples":
        ANALYSIS_BATCHES["muscle_F_101bp"]["samples"] +
        ANALYSIS_BATCHES["muscle_M_101bp"]["samples"] +
        ANALYSIS_BATCHES["muscle_F_76bp"]["samples"] +
        ANALYSIS_BATCHES["muscle_M_76bp"]["samples"],
    "sex": "both",
}

ANALYSIS_BATCHES["muscle_F"] = {
    "tissue": "muscle",
    "samples":
        ANALYSIS_BATCHES["muscle_F_101bp"]["samples"] +
        ANALYSIS_BATCHES["muscle_F_76bp"]["samples"],
    "sex": "F",
}

ANALYSIS_BATCHES["muscle_M"] = {
    "tissue": "muscle",
    "samples":
        ANALYSIS_BATCHES["muscle_M_101bp"]["samples"] +
        ANALYSIS_BATCHES["muscle_M_76bp"]["samples"],
    "sex": "M",
}

ANALYSIS_BATCHES["fibroblasts"] = {
    "tissue": "fibroblasts",
    "samples":
        ANALYSIS_BATCHES["fibroblasts_F"]["samples"] +
        ANALYSIS_BATCHES["fibroblasts_M"]["samples"],
    "sex": "both",
}
"""

