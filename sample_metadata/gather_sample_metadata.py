"""
This script collects RNA-seq sample metadata from several sources:
- Google sheets

- Terra workspaces used for running STAR and RNAseQC :
    https://app.terra.bio/#workspaces/macarthurlab-rnaseq-terra/macarthurlab-rnaseq-terra
    https://app.terra.bio/#workspaces/macarthurlab-rnaseq-terra/macarthurlab-rnaseq-terra%20-%20CMG_Broad_Orphan_Estonia-Onuap_RNA
-

"""


from __future__ import print_function
import os
import pandas as pd
import sys

import gspread

from google.oauth2.service_account import Credentials
from google.cloud import storage

from firecloud import api


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

#%%

# Get from google buckets

storage_client = storage.Client()
bucket = storage_client.get_bucket('macarthurlab-rnaseq')

#%%

sample_paths = {}

for p in bucket.list_blobs(): #prefix='*/'):
    print(p)

#%%
# Get from Google docs

creds = Credentials.from_service_account_file(
    '/Users/weisburd/project__rnaseq/code/rnaseq_methods/sample_metadata/seqr-project-fcfb2f45e8e4.json',
    scopes=['https://www.googleapis.com/auth/spreadsheets', "https://www.googleapis.com/auth/drive.file", "https://www.googleapis.com/auth/drive"]
)
client = gspread.authorize(creds)

#%%

sheets = client.open("RNA-seq metadata")

#%%
worksheets = sheets.worksheets()

#%%
all_samples_ws = sheets.get_worksheet(0)

#%%
beryls_ws = sheets.get_worksheet(1)

rows = beryls_ws.get()

beryls_ws_df = pd.DataFrame(data=rows[1:], columns=rows[0])


#%%

beryls_ws_df['Sample ID']

# Alias
# SpliceAI
# Viewer Link
# seqr link

#%%

import django
os.environ['DJANGO_SETTINGS_MODULE'] = 'settings'

django.setup()

#%%

"""
Not found:
38710
46N_RN
Dowling_7
Muntoni-1
Muntoni-2
"""
from seqr.models import Project, Family, Individual, Sample
print(Individual.objects.all().count())

sample_id_to_individual = {}
for sample_id in beryls_ws_df['Sample ID']:
    sample_id = sample_id.replace(".", "_")
    indivs = Individual.objects.filter(individual_id=sample_id)
    if not indivs:
        indivs = Individual.objects.filter(individual_id__contains=sample_id)
    if not indivs:
        indivs = Individual.objects.filter(individual_id=sample_id+"_1")
    if not indivs:
        indivs = Individual.objects.filter(individual_id=sample_id+"_M1")
    if len(indivs) > 1:
        print(sample_id, " ||| ".join(i.family.project.guid + " -- " + i.individual_id for i in indivs))

#sys.exit(0)

#%%
print(sheets.__dict__)
# Extract and print all of the values
list_of_hashes = sheet.get_all_records()
print(list_of_hashes)

sys.exit(0)


# RNA-seq metadata spreadsheet: https://docs.google.com/spreadsheets/d/1S3l28tZqFmzqqwqi_BCzuIkaVFmZz9eGpGtqtH5eVoo/edit#gid=363564316
RNASEQ_METADATA_SPREADSHEET_ID = '1S3l28tZqFmzqqwqi_BCzuIkaVFmZz9eGpGtqtH5eVoo'

def main():
    """Shows basic usage of the Sheets API.
    Prints values from a sample spreadsheet.
    """
    creds = None

    # The file token.pickle stores the user's access and refresh tokens, and is
    # created automatically when the authorization flow completes for the first
    # time.
    if os.path.exists('token.pickle'):
        with open('token.pickle', 'rb') as token:
            creds = pickle.load(token)

    # If there are no (valid) credentials available, let the user log in.
    if not creds or not creds.valid:
        if creds and creds.expired and creds.refresh_token:
            creds.refresh(Request())
        else:
            flow = InstalledAppFlow.from_client_secrets_file('credentials.json', SCOPES)
            creds = flow.run_local_server(port=0)

        # Save the credentials for the next run
        with open('token.pickle', 'wb') as token:
            pickle.dump(creds, token)

    service = build('sheets', 'v4', credentials=creds)

    # Call the Sheets API
    sheet = service.spreadsheets()
    sheet.va
    result = sheet.values().get(spreadsheetId=RNASEQ_METADATA_SPREADSHEET_ID).execute()
    values = result.get('values', [])

    if not values:
        print('No data found.')
    else:
        print('Name, Major:')
        for row in values:
            # Print columns A and E, which correspond to indices 0 and 4.
            print('%s, %s' % (row[0], row[4]))

if __name__ == '__main__':
    main()
