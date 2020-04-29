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
            os.path.expanduser('~/.config/gcloud/seqr-project-fcfb2f45e8e4.json'),
            scopes=[
                'https://www.googleapis.com/auth/spreadsheets',
                'https://www.googleapis.com/auth/drive.file',
                'https://www.googleapis.com/auth/drive',
            ]
        )

        _GSPREAD_CLIENT = gspread.authorize(creds)

    spreadsheet = _GSPREAD_CLIENT.open(spreadsheet_name)

    return spreadsheet


RNASEQ_METADATA_SPREADSHEET = get_spreasheet("RNA-seq metadata")
SEQR_INFO_AND_OTHER_METADATA_WORKSHEET = RNASEQ_METADATA_SPREADSHEET.worksheet("seqr info + other metadata (auto)")
DATA_PATHS_WORKSHEET = RNASEQ_METADATA_SPREADSHEET.worksheet("data paths (auto)")
BERYLS_WORKSHEET = RNASEQ_METADATA_SPREADSHEET.worksheet("Beryl's Supplementary Table 1")
BERYLS_WORKSHEET_2 = RNASEQ_METADATA_SPREADSHEET.worksheet("Beryl's RNAseq Probands")


def get_joined_metadata_df():
    df1_rows = DATA_PATHS_WORKSHEET.get()
    df1 = pd.DataFrame(data=df1_rows[1:], columns=df1_rows[0])
    df2_rows = SEQR_INFO_AND_OTHER_METADATA_WORKSHEET.get()
    df2 = pd.DataFrame(data=df2_rows[1:], columns=df2_rows[0])
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
