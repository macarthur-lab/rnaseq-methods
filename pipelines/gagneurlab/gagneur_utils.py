

GCLOUD_USER_ACCOUNT = "weisburd@broadinstitute.org"
GCLOUD_CREDENTIALS_LOCATION = "gs://weisburd-misc/creds"
GCLOUD_PROJECT = "seqr-project"

DOCKER_IMAGE = "weisburd/gagneurlab@sha256:be45788c8696a196bee25be269cb2de97277601bed65cbc3efcadc16acc5a764"
#DOCKER_IMAGE = "weisburd/gagneurlab@sha256:653d20def50adefa967587413eb7f8667272b67eeb90f935e3e372dda45dd367"

# https://i12g-gagneurweb.in.tum.de/public/workshops/RNAseq_ASHG19/input_data/annotations/gencode.v29lift37.annotation.txdb
GENCODE_TXDB = "gs://macarthurlab-rnaseq/gagneur/gencode.v26.annotation.txdb"
ALL_METADATA_TSV = "gs://macarthurlab-rnaseq/gagneur/metadata_table_for_all_RDG_and_GTEX_samples.tsv"
BAM_HEADER_PATH = "gs://macarthurlab-rnaseq/gagneur/fraser/bam_header.bam"

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
