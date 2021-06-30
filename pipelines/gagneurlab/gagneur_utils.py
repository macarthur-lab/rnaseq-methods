

GCLOUD_USER_ACCOUNT = "weisburd@broadinstitute.org"
GCLOUD_CREDENTIALS_LOCATION = "gs://weisburd-misc/creds"
GCLOUD_PROJECT = "seqr-project"

DOCKER_IMAGE = "weisburd/gagneurlab@sha256:be45788c8696a196bee25be269cb2de97277601bed65cbc3efcadc16acc5a764"
#DOCKER_IMAGE = "weisburd/gagneurlab@sha256:75a09c7ec42185b07206eb79f1d2f532ca712e19c18008779e9ea133153c807a"


# https://i12g-gagneurweb.in.tum.de/public/workshops/RNAseq_ASHG19/input_data/annotations/gencode.v29lift37.annotation.txdb
CLOUD_STORAGE_BASE_DIR = "gs://tgg-rnaseq/gagneur"
GENCODE_TXDB = f"{CLOUD_STORAGE_BASE_DIR}/gencode.v26.annotation.txdb"
#ALL_METADATA_TSV = f"{CLOUD_STORAGE_BASE_DIR}/metadata_table_for_all_RDG_and_GTEX_samples.tsv"
OUTRIDER_COUNTS_TSV_GZ = f"{CLOUD_STORAGE_BASE_DIR}/outrider/OUTRIDER_input_table_RDG_and_GTEX_counts_for_all_tissues.tsv.gz"
BAM_HEADER_PATH = f"{CLOUD_STORAGE_BASE_DIR}/fraser/bam_header.bam"

import gspread
import os
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


_OUTRIDER_RESULTS_SPREADSHEET = None
def get_OUTRIDER_results_spreadsheet():
    global _OUTRIDER_RESULTS_SPREADSHEET
    _OUTRIDER_RESULTS_SPREADSHEET = get_spreasheet("RNA-seq OUTRIDER results")
    return _OUTRIDER_RESULTS_SPREADSHEET


_FRASER_RESULTS_SPREADSHEET = None
def get_FRASER_results_spreadsheet():
    global _FRASER_RESULTS_SPREADSHEET
    _FRASER_RESULTS_SPREADSHEET = get_spreasheet("RNA-seq FRASER results")
    return _FRASER_RESULTS_SPREADSHEET


_RNASEQ_RESULTS_2020_12_10_SPREADSHEET = None
def get_RNASEQ_results_spreadsheet():
    global _RNASEQ_RESULTS_2020_12_10_SPREADSHEET
    _RNASEQ_RESULTS_2020_12_10_SPREADSHEET = get_spreasheet("RNA-seq results: FRASER, OUTRIDER")
    return _RNASEQ_RESULTS_2020_12_10_SPREADSHEET

_RNASEQ_TRUTH_DATA_SPREADSHEET = None
def get_rnaseq_truth_data_spreadsheet():
    global _RNASEQ_TRUTH_DATA_SPREADSHEET
    print("Loading 'RNA-seq truth data'")
    _RNASEQ_TRUTH_DATA_SPREADSHEET = get_spreasheet("RNA-seq truth data")
    return _RNASEQ_TRUTH_DATA_SPREADSHEET

