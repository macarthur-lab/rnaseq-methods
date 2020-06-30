
#%%

from __future__ import print_function
import collections
import json
import os
import pandas as pd
import re
import sys

from gspread_dataframe import set_with_dataframe

from sample_metadata.utils import get_seqr_info_and_other_metadata_worksheet, get_beryls_df, get_beryls_df_2, \
    get_data_paths_df, get_seqr_info_and_other_metadata_df

#%%
beryls_ws_df_merged = get_beryls_df().merge(get_beryls_df_2(), on="Sample ID", how="left")

"""
All columns in beryls_ws_df_merged:
Sample ID
Alias
Clinical Diagnosis
Sex
Age at muscle biopsy
Site of biopsy
Previous NGS testing
Status of genetic diagnosis before RNA-seq
Notes on genetic diagnosis status
Notes from paper
Data_type
Phenotype
Status
CanditateGenes (culprit,if solved)
Candidate  Variants
Variant  type
Variant consequence
Famliy data (all WES unless otherwise noted)
RNA BamPath
WGS BamPath
WES BAM
Data details
Collab PI
Realigned
Gender
Gender  confirmed
Tissue  verified
Sample  matched 
%Contamin  RNAseq
Ethnicity
Age at  Biopsy
Biopsy type
Phenotype comments
Other comments
Short Phenotype
Family
Include in manuscript?
"""

columns_from_beryls_worksheets = [
    "Sample ID",
    "Alias",
    "Clinical Diagnosis",
    "Sex",
    "Age at muscle biopsy",
    "Site of biopsy",
    #"Previous NGS testing",
    "Notes on genetic diagnosis status",
    "Notes from paper",
    "Data_type", # "sample types (Beryl)",
    "Phenotype", # "phenotype (Beryl)",
    "Status",  # "solved? (Beryl)",
    'CanditateGenes\n(culprit,if solved)',
    'Candidate \nVariants',
    'Variant \ntype',
    'Variant consequence',
    'Data details',
    'Collab PI',
    '%Contamin\n RNAseq',
    'Ethnicity',
    'Age at \nBiopsy',
    'Biopsy type',
    'Phenotype comments',
    'Other comments',
    'Short Phenotype',
    'Include in manuscript?',
]
beryls_ws_df_merged = beryls_ws_df_merged[columns_from_beryls_worksheets]

beryls_ws_df_merged = beryls_ws_df_merged.rename(columns={
    "Status": "Genetic diagnosis Status",
})

beryls_ws_df_merged = beryls_ws_df_merged.rename(columns={c: c.replace("\n", " ") + " (Beryl)" for c in columns_from_beryls_worksheets if c != "Sample ID"})

beryls_ws_df_merged["Notes (Beryl)"] = beryls_ws_df_merged["Notes on genetic diagnosis status (Beryl)"].fillna('')

for i, row in beryls_ws_df_merged.iterrows():
    #print(i, ",", df_export.at[i, "Notes (Beryl)"], ",", row["Notes from paper"])
    if row["Notes from paper (Beryl)"] and not isinstance(row["Notes from paper (Beryl)"], float):
        beryls_ws_df_merged.at[i, "Notes (Beryl)"] += "\n" + row["Notes from paper (Beryl)"]

# %%

# join DATA_PATHS_WORKSHEET with BERYLS_WORKSHEET

data_paths_ws_df = get_data_paths_df()
beryls_ws_df_merged['sample_id'] = ""
for s in beryls_ws_df_merged['Sample ID']:
    all_samples_row = data_paths_ws_df().loc[data_paths_ws_df['sample_id'].str.contains(s), ]
    if all_samples_row.shape[0] == 0:
        print("sample id " + s + " from Beryl's table not found in data_paths_ws_df")
        continue
    elif all_samples_row.shape[0] > 1:
        print("sample id " + s + " from Beryl's table matched more than 1 entry in data_paths_ws_df: " + ", ".join(all_samples_row.sample_id.tolist()))
        continue

    beryls_ws_df_merged.loc[beryls_ws_df_merged['Sample ID'] == s, 'sample_id'] = all_samples_row.sample_id.item()

joined_df = data_paths_ws_df.merge(beryls_ws_df_merged, on="sample_id", how="left")
joined_df2 = beryls_ws_df_merged.merge(data_paths_ws_df, on="sample_id", how="left")

#joined_df
print("Found match for", joined_df2[~joined_df2['hg19_bam'].isnull()].shape[0], "out of", beryls_ws_df_merged.shape[0], "rows")
print("Match not found for\n  '" + "'\n  '".join(joined_df2[joined_df2['hg19_bam'].isnull()]['Sample ID']))

print("Found match for", joined_df[~joined_df['Sex (Beryl)'].isnull()].shape[0], "out of", beryls_ws_df_merged.shape[0], "rows")
print("-----------------")
print("joined_df columns:")
print('  "'+ '",\n  "'.join(joined_df.columns) + '",')


#%%

# update batch_date_from_hg19_bam_header column

from sample_metadata.utils import get_date_from_bam_header

seqr_info_and_other_metadata_ws_rows = get_seqr_info_and_other_metadata_df()
joined_df3 = joined_df.merge(seqr_info_and_other_metadata_ws_rows[["sample_id", "batch_date_from_hg19_bam_header"]], on="sample_id", how="left")
for i, row in joined_df3.iterrows():
    if not row["batch_date_from_hg19_bam_header"]:
        d = joined_df3.at[i, "batch_date_from_hg19_bam_header"] = get_date_from_bam_header(row["hg19_bam"])
        print(d, joined_df3.at[i, "hg19_bam"])

joined_df3

#%%

# parse rnaseqc metrics.txt

from utils import get_rnaseqc_metrics

RNASEQC_COLUMNS =  ["stranded? (rnaseqc)", "read length (rnaseqc)", "total reads x 10^6 (rnaseqc)", "mapping rate (rnaseqc)"]
joined_df4 = joined_df3.merge(seqr_info_and_other_metadata_ws_rows[["sample_id"] + RNASEQC_COLUMNS], on="sample_id", how="left")
for i, row in joined_df4.iterrows():
    rnaseqc_metrics_file_path = row["rnaseqc_metrics"].strip()
    if not rnaseqc_metrics_file_path or row["stranded? (rnaseqc)"]:
        continue

    #print(row["rnaseqc_metrics"])
    metrics_dict = get_rnaseqc_metrics(rnaseqc_metrics_file_path)
    end1_sense_rate = float(metrics_dict['End 1 Sense Rate'])
    end2_sense_rate = float(metrics_dict['End 2 Sense Rate'])

    if 0.45 < end1_sense_rate < 0.55 and 0.45 < end2_sense_rate < 0.55:
        stranded = "no"
    elif end1_sense_rate < 0.05 and end2_sense_rate > 0.95:
        stranded = "yes"
    else:
        stranded = None
        print("ERROR: End 1 Sense Rate and/or End 2 Sense Rate are out of bounds: %0.2f %0.2f" % (end1_sense_rate, end2_sense_rate))

    print("%0.2f   %0.2f" % (float(metrics_dict['End 1 Sense Rate']), float(metrics_dict['End 2 Sense Rate'])))
    joined_df4.at[i, "stranded? (rnaseqc)"] = stranded
    joined_df4.at[i, "read length (rnaseqc)"] = "%0.0f" % float(metrics_dict['Read Length'])
    joined_df4.at[i, "total reads x 10^6 (rnaseqc)"] = "%0.0f" % (int(metrics_dict['Total Reads'])/10**6)
    joined_df4.at[i, "mapping rate (rnaseqc)"] = "%0.2f" % float(metrics_dict['Mapping Rate'])

    #print("stranded?", joined_df4.at[i, "stranded? (rnaseqc)"], rnaseqc_metrics_file_path)
    print("read length", joined_df4.at[i, "read length (rnaseqc)"], rnaseqc_metrics_file_path)

joined_df4


# Alias
# SpliceAI
# Viewer Link
# seqr link

#%%

sys.path.append(os.path.expanduser("~/code/seqr"))

import django
os.environ['DJANGO_SETTINGS_MODULE'] = 'settings'

django.setup()


#%%

"""
            if i.phenotips_data:
                phenotips_json = json.loads(i.phenotips_data)
                phenotips_fields = _parse_phenotips_data(phenotips_json)
            else:
                phenotips_fields = {}

            if include_hpo_terms_present:
                row.append(phenotips_fields.get('phenotips_features_present', ''))
            if include_hpo_terms_absent:
                row.append(phenotips_fields.get('phenotips_features_absent', ''))
            if include_paternal_ancestry:
                row.append(phenotips_fields.get('paternal_ancestry', ''))
            if include_maternal_ancestry:
                row.append(phenotips_fields.get('maternal_ancestry', ''))
            if include_age_of_onset:
                row.append(phenotips_fields.get('age_of_onset', ''))

"""
"""
Not found:
46N_RN
Dowling_7
Muntoni-1
Muntoni-2
"""

sample_id_to_seqr_indiv_id = {

    #"SWE_36223_1", 	 # control
    #"SWE_36326_1",      # following up on potential novel gene
    #"SRPK3_EM__1",
    #"SRPK3_JAP_1",
    #"SRPK3_USA_1",

    "431-1_A_3":  "GAZ_431_3_D2",
    "431-2_A_2":  "GAZ_431_2_D1",  # unaffected mother?
    "432-3_A_2":  "GAZ_432_3_D1",
    "437-2_A_2":  "GAZ_437_2_D1",  # unaffected mother?
    "437-3_A_2":  "GAZ_437_3_D1",
    "447-3_A_3":  "GAZ_447_3_D2",

    "BEG_1435-01_T1370": "BEG_1435-1_01",
    "BEG_1111-1_T1105": "BEG_1111-1_01",
    "BEG_1025-1_T999": "BEG_1025-1_1",
    "BEG_1078-1_T1071": "BEG_1078-1_01",
    "BEG_1230-1_T1227": "BEG_1230-1",
    "BEG_14-4_T65": "BEG_14-4_1",
    "BEG_1438-1_T1339": "BEG_1438-1_01",
    "BEG_851-1_T840": "BEG_851-1_1",
    "BEG_887-1_T1240": "BEG_887-1_1",
    "BEG_916-1_T916": "BEG_916-1_1",

    "BON_B09-27-1_1":   "B09-27-1",
    "BON_B12-33-2_1":   "B12-33-2",
    "BON_B12-74-1_1":   "B12-74-1",
    "BON_B12-76-1_2":   "B12-76-1",
    "BON_B13-55_1_2":   "BON_B13-55_1_1",
    "BON_B14-163-1_2":  "B14-163-1",
    "BON_B14-20_1":     "B14-20",
    "BON_B14-51_1_2":   "BON_B14-51_1_1",
    "BON_B14-60-1_2":   "B14-60-1",
    "BON_B14-71-2_1":   "B14-71-2",
    "BON_B14-75-1_1":   "B14-75-1",
    "BON_B15-118_1":    "BON_B15-118_2",
    "BON_B15-125_1_2":  "BON_B15-125_1_1",
    "BON_B15-26_1_1":   "BON_B15-26_1",
    "BON_B15-76-2_2":   "B15-76-2",
    "BON_B15-98_1_2":   "BON_B15-98_1_1",
    "BON_B16-19_1":     "BON_B16-19_2",
    "BON_B16-22_1":     "BON_B16-22-1",
    "BON_B16-50_1_2":   "BON_B16-50_1_1",
    "BON_B16-53-1_1":   "BON_B16-53-1",
    "BON_B16-57_1_2":   "BON_B16-57_1_1",
    "BON_B16-75-1_2":   "BON_B16-75_1_1",
    "BON_UC219-1_1":    "UC219-1",

    "HK-071-001": "HK071-001",
    "HK-071-004": "HK071-004",
    "HK-006": "HK006_0016",
    "HK-010": "HK010_0026",
    "HK-015": "HK015_0036",
    "HK-022": "HK022_0059",
    "HK-025": "HK025_0068",
    "HK-032": "HK032_0081",
    "HK-035": "HK035_0088",
    "HK-047": "HK047_0116",
    "HK-061": "HK061-0157",
    "HK-070": "HK070-0180",
    "HK-073": "HK073-001",
    "HK-075": "HK075-001",
    "HK-081": "HK081-001",
    "HK-085": "HK085-001",  ### ???
    "HK-087": "HK087-001",
    "HK-095": "HK095-001_1",
    "HK-100": "HK100-001",
    "HK-104": "HK104-001_D2",
    "HK-107": "HK107-001",
    "HK-108": "HK108-001",
    "HK-114": "HK114-001",
    "HK-115": "HK115-001",
    "HK-116": "HK116-001",
    "HK-119": "HK119-001",
    "HK-124": "OUN_HK124_001",
    "HK-126": "OUN_HK126_001",
    "HK-127": "OUN_HK127-001_D1",
    "HK-131": "OUN_HK131_001",
    "HK-133": "OUN_HK133-001_D1",
    "HK-134": "OUN_HK134_001",
    "HK-139": "OUN_HK139_001",

    "HK018_0047_2": "HK018_0047",
    "HK069-0177_2": "HK069-0177",
    "HK072-001_2": "HK072-001",
    "HK088-001_2": "HK088-001",


    "RGP_554_3_2": "RGP_554_3",

    "RGP_54_3_2": "RGP_54_3",
    "RGP_56_3_3": "RGP_56_3",
    "RGP_7_1_2": "RGP_7_1",
    "RGP_7_2_2": "RGP_7_2",
    "RGP_7_3_2": "RGP_7_3",
    "RGP_7_4_2": "RGP_7_4",
    "RGP_7_5_2": "RGP_7_5",
    "RGP_800_3_2": "RGP_800_3",
    "RGP_94_3_2": "RGP_94_3",
    "RGP_94_4_2": "RGP_94_4",

    "OUN_HK018_0047": "HK018_0047",
    "OUN_HK047_0116": "HK047_0116",
    "OUN_HK079_001": "HK079-001",
    "OUN_HK080_001": "HK080-001",
    "OUN_HK081_001": "HK081-001",
    "OUN_HK112_001": "HK112-001",
    "OUN_HK116_001": "HK116-001",
    "OUN_HK137_001": "OUN_HK137-001_D1",

    "MAN_1438-01-M1": "MAN_1438-01_1",
    "MAN_1001_01_M1_D1": "MAN_1001-01_1",
    "ICCV_458_10CC06258_02": "ICCV_458_10CC06258_01",

    "VCGS_DLS_FAM1_1": "VCGS_FAM1_1",
    "VCGS_FAM11_32_2": "VCGS_FAM11_32",
    "VCGS_FAM12_36_2": "VCGS_FAM12_36",
    "VCGS_FAM147_459_2": "VCGS_FAM147_459",
    "VCGS_FAM148_462_2": "VCGS_FAM148_462",
    "VCGS_FAM149_465_2": "VCGS_FAM149_465",
    "VCGS_FAM150_468_2": "VCGS_FAM150_468",
    "VCGS_FAM1_1_2": "VCGS_FAM1_1",
    "VCGS_FAM22_69_2": "VCGS_FAM22_69",
    "VCGS_FAM26_81_2": "VCGS_FAM26_81",
    "VCGS_FAM27_84_2": "VCGS_FAM27_84",
    "VCGS_FAM2_4_2": "VCGS_FAM2_4",
    "VCGS_FAM31_97_2": "VCGS_FAM31_97",
    "VCGS_FAM3_7_2": "VCGS_FAM3_7",
    "VCGS_FAM42_132_2": "VCGS_FAM42_132",
    "VCGS_FAM4_13_2": "VCGS_FAM4_13",
    "VCGS_FAM52_162_3": "VCGS_FAM52_162",
    "VCGS_FAM73_227_2": "VCGS_FAM73_227",

    "B13-07-1RNA": "B13-07-1",
    "B15-25_1_2": "B15-25_1",
    "B15-28_1_1": "B15-28_1",
    "B16-47_1_2": "B16-47_1",
    "BB0280_CH_AffF_2": "BB0290-CH-AffC",

    "NH12-843_Fibroblasts": "NH12-843",
    "NH12-843_MyoD_Day5": "NH12-843",

    "TOP_MAAC031_F_Muscle1": "MAAC031",
    "MBEL028_002_1": "MBEL028_2",

    "MBRU030_2": "MBRU030",
    "MESP014_2": "MESP014",
    "MESP039_2": "MESP039",
    "MESP021_001_2": "MESP021_001",
    "MBEL028_001_3": "MBEL028_001",
    "MCOP008_001_2": "MCOP008_001",
    "MGLA003_001_2": "MGLA003_001",
    "MMAD002_001_2": "MMAD002_001",
    "MTEH041_001_2": "MTEH041_001",

    "MAN_0063_01_02": "MAN_0063_01_01",  # deceased - 2 different tissues from same autopsy, degrade low quality RNA
    "MAN_0063_01_03": "MAN_0063_01_01",  # deceased - 2 different tissues from same autopsy, degrade low quality RNA

    "LIA_TIS03_2": "LIA_TIS03_1",
    "LIA_MAS02_2": "LIA_MAS02",
    "K1157-1-4": "K1157-1",
}


no_seqr_record = {
    "SWE_36223_1",
    "SWE_36326_1",
    "BON_B18-25_1",
    "SRPK3_EM__1",
    "SRPK3_JAP_1",
    "SRPK3_USA_1",
    "431-1_A_3",
    "431-2_A_2",
    "432-3_A_2",
    "437-2_A_2",
    "437-3_A_2",
    "447-3_A_3",
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
    "MUN_FAM1_CTRL_01",
    "MUN_FAM1_CTRL_02",
    "MUN_FAM1_CTRL_03",
    "MUN_FAM1_CTRL_04",
    "MUN_FAM2_TOTALMDC1A1_02",
    "MUN_FAM2_TOTALMDC1A1_03",
    "MUN_FAM2_TOTALMDC1A1_04",
    "MUN_FAM2_TOTALMDC1A1_05",
    "MUN_FAM3_PARTIALMDC1A1_01",
    "MUN_FAM3_PARTIALMDC1A1_02",
    "MUN_FAM3_PARTIALMDC1A1_03",
    "MUN_FAM3_PARTIALMDC1A1_05",
    "MUN_FAM4_ATYPICALMDC1A_01_R1",
    "MUN_FAM4_ATYPICALMDC1A_01_R2",
    "MUN_FAM4_ATYPICALMDC1A_01_R3",
    "MUN_FAM5_SIBLINGMDC1A_01_R1",
    "MUN_FAM5_SIBLINGMDC1A_01_R2",
    "MUN_FAM5_SIBLINGMDC1A_01_R3",
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


from seqr.models import Project, Family, Individual, Sample, IgvSample, SavedVariant, VariantTag, VariantNote
#print(Individual.objects.all().count())

counters = collections.defaultdict(int)
sample_id_to_indivs = {}
for sample_id in joined_df4['sample_id']:
    if sample_id in no_seqr_record:
        continue

    # find matching indivs
    seqr_id = sample_id_to_seqr_indiv_id.get(sample_id, sample_id)
    indivs = Individual.objects.filter(individual_id=seqr_id)
    if not indivs:
        indivs = Individual.objects.filter(individual_id__contains=seqr_id)
    if not indivs:
        indivs = Individual.objects.filter(individual_id=seqr_id+"_1")
    if not indivs:
        indivs = Individual.objects.filter(individual_id=seqr_id+"_M1")
    if len(indivs) > 1:
        pass
        #print(sample_id + ":", " ||| ".join(i.family.project.name + "/" + i.family.family_id for i in indivs))
    elif len(indivs) == 0:
        print(sample_id, type(sample_id))

    counters[len(indivs)] += 1

    assert len(indivs) <= 2, "Expected no more than 2 records per sample_id"

    sample_id_to_indivs[sample_id] = indivs

print("---")
print(counters)

#%%

analysis_status_lookup = dict(Family.ANALYSIS_STATUS_CHOICES)

seqr_fields_by_sample_id = collections.defaultdict(lambda: collections.defaultdict(list))
for sample_id, indivs in sample_id_to_indivs.items():

    seqr_fields = seqr_fields_by_sample_id[sample_id]
    seqr_fields['sample_id'] = sample_id
    variant_tag_counters = collections.defaultdict(int)
    num_variant_notes = 0
    for indiv in indivs:
        family = indiv.family
        project = family.project

        saved_variant_ids_for_family = set([sv.pk for sv in family.savedvariant_set.all()])

        for k, n in collections.Counter([tag.varianttag.variant_tag_type.name for tag in VariantTag.saved_variants.through.objects.filter(savedvariant_id__in=saved_variant_ids_for_family)]).items():
            variant_tag_counters[k] += n

        num_variant_notes += len([note for note in VariantNote.saved_variants.through.objects.filter(savedvariant_id__in=saved_variant_ids_for_family)])

        #[sv.varianttag_set.all() for sv in family.savedvariant_set.all()]
        #[sv.variantnote_set.all() for sv in family.savedvariant_set.all()]

        if not seqr_fields['proj (seqr)']:
            project_i = ""
        elif not seqr_fields['proj2 (seqr)']:
            project_i = 2
        else:
            print("ERROR: more than 2 projects have individual: " + str(indiv) + ". Skipping...")
            continue

        project_page_url = "https://seqr.broadinstitute.org/project/%s/project_page" % (project.guid)
        family_page_url = "https://seqr.broadinstitute.org/project/%s/family_page/%s" % (project.guid, family.guid)

        seqr_fields['proj%s (seqr)' % project_i] = '=HYPERLINK("%s", "%s")' % (project_page_url, project.name)
        seqr_fields['fam%s (seqr)' % project_i] = '=HYPERLINK("%s", "%s")' % (family_page_url, family.family_id)

        seqr_fields['genome (seqr)'].append("hg%s" % project.genome_version)

        if family.analysis_status:
            seqr_fields['analysis status (seqr)'].append(analysis_status_lookup[family.analysis_status])

        if family.coded_phenotype:
            seqr_fields['coded phenotype (seqr)'].append(family.coded_phenotype)
        if family.analysis_summary or family.analysis_notes:
            seqr_fields['anlaysis summary + notes (seqr)'].append((family.analysis_summary or "" + "\n" + family.analysis_notes or "").strip())
        if family.internal_case_review_notes:
            seqr_fields['internal case review notes (seqr)'].append((family.internal_case_review_notes or "").strip())

        if indiv.sex and indiv.sex != "U":
            seqr_fields['sex'].append(indiv.sex)
        if indiv.population:
            seqr_fields['population (seqr)'].append(indiv.population or "")

        sample_types = [sample.sample_type for sample in indiv.sample_set.all() if sample.sample_type and sample.sample_type.strip() and sample.sample_type != "RNA"]
        if sample_types:
            seqr_fields['sample type (seqr)'].append((" ".join(set(sample_types))))
        cram_paths = [sample.file_path for sample in indiv.igvsample_set.all()]
        if cram_paths:
            seqr_fields['cram path (seqr)'].append(" ".join(set(cram_paths)))

        #if indiv.phenotips_data:
        #    print(json.loads(indiv.phenotips_data))

        # indiv.display_name
        #if indiv.mother:
        #    indiv.mother.individual_id
        #if indiv.father:
        #    indiv.father.individual_id
        # sample.dataset_type
        # sample.elasticsaerch_index


    if joined_df4[joined_df4['sample_id'] == sample_id].shape[0] > 0:
        other_sex_annotation = joined_df4[joined_df4['sample_id'] == sample_id]['Sex (Beryl)']
        if type(other_sex_annotation) != float and type(other_sex_annotation.item()) != float and other_sex_annotation.item():
            seqr_fields['sex'].append(other_sex_annotation.item())

    seqr_fields['sex'] = ", ".join(sorted(set(seqr_fields['sex']))).encode("utf-8")
    seqr_fields['genome (seqr)'] = ", ".join(sorted(set(seqr_fields['genome (seqr)']))).encode("utf-8")
    seqr_fields['population (seqr)'] = ", ".join(sorted(set(seqr_fields['population (seqr)']))).encode("utf-8")
    seqr_fields['sample type (seqr)'] = ", ".join(sorted(set(seqr_fields['sample type (seqr)']))).encode("utf-8")
    seqr_fields['analysis status (seqr)'] = ", ".join(sorted(set(seqr_fields['analysis status (seqr)']))).encode("utf-8")

    seqr_fields['variant tags (seqr)'] = ", ".join(["%s %s tag(s)" % (n, key) for key, n in sorted(variant_tag_counters.items(), key=lambda x: x[0], reverse=True)]) + " "
    seqr_fields['variant notes (seqr)'] = "%s" % num_variant_notes

    seqr_fields['coded phenotype (seqr)'] = ", ".join(seqr_fields['coded phenotype (seqr)']).encode("utf-8") + " "
    seqr_fields['anlaysis summary + notes (seqr)'] = "\n".join(seqr_fields['anlaysis summary + notes (seqr)']).encode("utf-8") + " "
    seqr_fields['internal case review notes (seqr)'] = "\n".join(seqr_fields['internal case review notes (seqr)']).encode("utf-8") + " "
    seqr_fields['cram path (seqr)'] = ", ".join(sorted(set(seqr_fields['cram path (seqr)']))).encode("utf-8") + " "

seqr_fields_by_sample_id.values()

#sys.exit(0)
#%%
SEQR_INFO_COLUMNS = [
    'proj (seqr)',
    'fam (seqr)',
    'proj2 (seqr)',
    'fam2 (seqr)',
    'sex',
    'genome (seqr)',
    'population (seqr)',
    'sample type (seqr)',
    'analysis status (seqr)',
    'variant tags (seqr)',
    'variant notes (seqr)',
    'coded phenotype (seqr)',
    'anlaysis summary + notes (seqr)',
    'internal case review notes (seqr)',
    'cram path (seqr)',
]

seqr_info_df = pd.DataFrame(
    columns=["sample_id"] + SEQR_INFO_COLUMNS,
    data=seqr_fields_by_sample_id.values()
)

seqr_info_df


#%%

joined_df5 = joined_df4.merge(seqr_info_df, on="sample_id", how="left")

joined_df5["solved using RNA-seq?"] = ""
joined_df5["solved not using RNA-seq?"] = ""

# %%

df_export = joined_df5[[
    "sample_id",
    "star_pipeline_batch",
    "batch_date_from_hg19_bam_header",
] + RNASEQC_COLUMNS + [
    "solved using RNA-seq?",
    "solved not using RNA-seq?"
] + SEQR_INFO_COLUMNS + [c for c in beryls_ws_df_merged.columns if c not in ('sample_id', 'Sample ID')]]

# export joined data to SEQR_INFO_AND_OTHER_METADATA_WORKSHEET
set_with_dataframe(SEQR_INFO_AND_OTHER_METADATA_WORKSHEET, df_export.fillna(''), resize=True)

print("Updated", SEQR_INFO_AND_OTHER_METADATA_WORKSHEET.title)


# %%
