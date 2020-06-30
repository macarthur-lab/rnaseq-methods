import gspread
import os
import pandas as pd
import subprocess
from google.oauth2.service_account import Credentials


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
_RNASEQ_METADATA_SPREADSHEET = None
_SEQR_INFO_AND_OTHER_METADATA_WORKSHEET = None
_DATA_PATHS_WORKSHEET = None
_IMPUTED_METADATA_WORKSHEET = None
_BERYLS_WORKSHEET = None
_BERYLS_WORKSHEET_2 = None

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


def get_seqr_info_and_other_metadata_worksheet():
    global _SEQR_INFO_AND_OTHER_METADATA_WORKSHEET
    _SEQR_INFO_AND_OTHER_METADATA_WORKSHEET = get_rnaseq_metadata_spreadsheet().worksheet("seqr info + other metadata (auto)")
    return _SEQR_INFO_AND_OTHER_METADATA_WORKSHEET


def get_data_paths_worksheet():
    global _DATA_PATHS_WORKSHEET
    _DATA_PATHS_WORKSHEET = get_rnaseq_metadata_spreadsheet().worksheet("data paths (auto)")
    return _DATA_PATHS_WORKSHEET


def get_imputed_metadata_worksheet():
    global _IMPUTED_METADATA_WORKSHEET
    _IMPUTED_METADATA_WORKSHEET = get_rnaseq_metadata_spreadsheet().worksheet("imputed (auto)")
    return _IMPUTED_METADATA_WORKSHEET


def get_beryls_worksheet():
    global _BERYLS_WORKSHEET
    _BERYLS_WORKSHEET = get_rnaseq_metadata_spreadsheet().worksheet("Beryl's Supplementary Table 1")
    return _BERYLS_WORKSHEET


def get_beryls_worksheet_2():
    global _BERYLS_WORKSHEET_2
    _BERYLS_WORKSHEET_2 = get_rnaseq_metadata_spreadsheet().worksheet("Beryl's RNAseq Probands")
    return _BERYLS_WORKSHEET_2


def get_seqr_info_and_other_metadata_df():
    rows = get_seqr_info_and_other_metadata_worksheet().get()
    return pd.DataFrame(data=rows[1:], columns=rows[0])


def get_data_paths_df():
    rows = get_data_paths_worksheet().get()
    return pd.DataFrame(data=rows[1:], columns=rows[0])


def get_imputed_metadata_df():
    rows = get_imputed_metadata_worksheet().get()
    return pd.DataFrame(data=rows[1:], columns=rows[0])


def get_beryls_df():
    rows = get_beryls_worksheet().get()
    return pd.DataFrame(data=rows[1:], columns=rows[0])


def get_beryls_df_2():
    rows = get_beryls_worksheet_2().get()
    return pd.DataFrame(data=rows[1:], columns=rows[0])


def get_joined_metadata_df():
    df1 = get_data_paths_df()
    df2 = get_seqr_info_and_other_metadata_df()
    df2 = df2[[c for c in df2.columns if c not in ("star_pipeline_batch")]]  # remove columns that exist in both tables
    return df1.merge(df2, on="sample_id", how="left").set_index("sample_id", drop=False)


def get_date_from_bam_header(bam_path):
    output = subprocess.check_output("gsutil cat %s | samtools view -H - | grep ^@RG | head -n 1" % bam_path, shell=True)
    read_group_annotations = {}
    for i, rg_field in enumerate(output.rstrip().split("\t")):
        if i == 0:
            continue  # skip @RG prefix
        key = rg_field[:rg_field.find(":")]
        value = rg_field[rg_field.find(":")+1:]
        read_group_annotations[key] = value

    bam_date = read_group_annotations['DT'][:7]
    return bam_date


def get_rnaseqc_metrics(rnaseqc_metrics_file_path):
    output = subprocess.check_output("gsutil cat %s" % rnaseqc_metrics_file_path, shell=True)
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


def get_gtex_rnaseq_sample_metadata_df():
    rows = get_gtex_rnaseq_sample_metadata_worksheet().get()
    return pd.DataFrame(data=rows[1:], columns=rows[0]).set_index("SAMPID", drop=False)


def get_gtex_wes_sample_metadata_df():
    rows = get_gtex_wes_sample_metadata_worksheet().get()
    return pd.DataFrame(data=rows[1:], columns=rows[0]).set_index("SAMPID", drop=False)


def get_gtex_wgs_sample_metadata_df():
    rows = get_gtex_wgs_sample_metadata_worksheet().get()
    return pd.DataFrame(data=rows[1:], columns=rows[0]).set_index("SAMPID", drop=False)
