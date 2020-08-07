
#%%

from __future__ import print_function
import collections
import os
import pandas as pd
import sys

pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)

from gspread_dataframe import set_with_dataframe

from sample_metadata.utils import \
    get_seqr_info_and_other_metadata_worksheet, \
    get_beryls_supplementary_table_df, \
    get_beryls_rnaseq_probands_df, \
    get_beryls_seqr_data_df, \
    get_data_paths_df, \
    get_seqr_info_and_other_metadata_df, \
    get_rnaseqc_metrics, \
    get_date_from_bam_header


#%%

def left_join_beryls_table(df, beryls_df):
    """Joins df 'sample_id' column to beryls_df 'Sample ID' column"""
    beryls_df = beryls_df.copy()
    beryls_df['sample_id'] = ""

    found_sample_ids = set()
    for s in beryls_df['Sample ID']:
        all_samples_row = df.loc[df['sample_id'].str.contains(s), ]
        if all_samples_row.shape[0] == 0:
            print("sample id " + s + " from Beryl's table not found in data_paths_df")
            continue
        elif all_samples_row.shape[0] > 1:
            print("sample id " + s + " from Beryl's table matched more than 1 entry in data_paths_ws_df: " + ", ".join(all_samples_row.sample_id.tolist()))
            continue

        beryls_df.loc[beryls_df['Sample ID'] == s, 'sample_id'] = all_samples_row.sample_id.item()
        found_sample_ids.add(s)

    not_found_sample_ids = set(beryls_df['Sample ID']) - found_sample_ids
    print(f"{len(found_sample_ids)} out of {len(beryls_df)} sample ids found in left side df: {not_found_sample_ids}")
    print(f"{len(not_found_sample_ids)} sample ids from right-side df not found in left side df: {not_found_sample_ids}")

    joined_df = df.merge(beryls_df, on="sample_id", how="left")

    print("-----------------")
    print("joined_df columns:")
    print('  "' + '",  "'.join(joined_df.columns) + '",')

    return joined_df


#%%

data_paths_df = get_data_paths_df()
print(data_paths_df.shape)
print(data_paths_df.columns)

print("Number of non-null data paths: ")
get_count = lambda c: sum(data_paths_df[c].str.len() > 0)
print("\n".join([f"{get_count(c)} {c}" for c in sorted(data_paths_df.columns, key=get_count, reverse=True)]))


#%%

final_df = data_paths_df[['sample_id',	'star_pipeline_batch', "hg19_bam", "rnaseqc_metrics"]]

#%%
# download current table to allow partial updates to some columns
seqr_info_and_other_metadata_rows = get_seqr_info_and_other_metadata_df()

final_df = final_df.merge(seqr_info_and_other_metadata_rows[[
    "sample_id",
    "batch_date_from_hg19_bam_header",
    "imputed tissue",
    "imputed sex"
]], on="sample_id", how="left")


#%%
# update batch_date_from_hg19_bam_header column

for i, row in final_df.iterrows():
    if not row.batch_date_from_hg19_bam_header or isinstance(row.batch_date_from_hg19_bam_header, float):
        print(row.hg19_bam)
        d = get_date_from_bam_header(row.hg19_bam)
        print(d)
        final_df.at[i, "batch_date_from_hg19_bam_header"] = d

#%%

# parse rnaseqc metrics

RNASEQC_COLUMNS = [
    "stranded? (rnaseqc)",
    "read length (rnaseqc)",
    "total reads x 10^6 (rnaseqc)",
    "mapping rate (rnaseqc)"
]

final_df = final_df.merge(seqr_info_and_other_metadata_rows[["sample_id"] + RNASEQC_COLUMNS], on="sample_id", how="left")
for i, row in final_df.iterrows():
    rnaseqc_metrics_file_path = row["rnaseqc_metrics"].strip()
    if not rnaseqc_metrics_file_path:
        continue

    if row["stranded? (rnaseqc)"] and not isinstance(row["stranded? (rnaseqc)"], float):
        # skip rows that already have rnaseqc values
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
    final_df.at[i, "stranded? (rnaseqc)"] = stranded
    final_df.at[i, "read length (rnaseqc)"] = "%0.0f" % float(metrics_dict['Read Length'])
    final_df.at[i, "total reads x 10^6 (rnaseqc)"] = "%0.0f" % (int(metrics_dict['Total Reads'])/10**6)
    final_df.at[i, "mapping rate (rnaseqc)"] = "%0.2f" % float(metrics_dict['Mapping Rate'])

    #print("stranded?", final_df.at[i, "stranded? (rnaseqc)"], rnaseqc_metrics_file_path)
    print("read length", final_df.at[i, "read length (rnaseqc)"], rnaseqc_metrics_file_path)

final_df


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

    "WAL_OTH2405f-Dec19",
    "WAL_OTH2405f-Feb20A",
    "WAL_OTH2405f-Feb20C",
    "WAL_OTH2406d-Dec19",
    "WAL_OTH2406d-Feb20A",
    "WAL_OTH2406d-Feb20B",
    "WAL_OTH2407d-Dec19",
    "WAL_OTH2407d-Feb20A",
    "WAL_OTH2407d-Feb20B",
    "WAL_OTH2419f-Dec19",
    "WAL_OTH2419f-Feb20A",
    "WAL_OTH2419f-Feb20B",
}


#%%
# populate sample_id_to_indivs dict

from seqr.models import Project, Family, Individual, Sample, IgvSample, SavedVariant, VariantTag, VariantNote
#print(Individual.objects.all().count())

sample_id_to_indivs = {}

counters = collections.defaultdict(int)
for sample_id in final_df['sample_id']:
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
        print(f"seqr indiv not found for: {sample_id}")

    counters[len(indivs)] += 1

    assert len(indivs) <= 2, "Expected no more than 2 records per sample_id"

    sample_id_to_indivs[sample_id] = indivs

print("---")
print(counters)

#%%

# add seqr metadata columns

analysis_status_lookup = dict(Family.ANALYSIS_STATUS_CHOICES)

seqr_fields_by_sample_id = collections.defaultdict(lambda: collections.defaultdict(list))
for sample_id, indivs in sample_id_to_indivs.items():

    seqr_fields = seqr_fields_by_sample_id[sample_id]
    seqr_fields['sample_id'] = sample_id
    seqr_fields['indiv (seqr)'] = indivs[0].individual_id
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
            seqr_fields['sex'].append(indiv.sex.decode('UTF-8') if isinstance(indiv.sex, bytes) else indiv.sex)
        if indiv.population:
            seqr_fields['population (seqr)'].append(indiv.population.decode('UTF-8') if isinstance(indiv.population, bytes) else indiv.population)

        samples = [sample for sample in indiv.sample_set.all() if sample.sample_type and sample.sample_type.strip() and sample.sample_type != "RNA"]

        sample_ids = [(sample.sample_id.decode('UTF-8') if isinstance(sample.sample_id, bytes) else sample.sample_id) for sample in samples]
        if sample_ids:
            seqr_fields['sample id (seqr)'].append((", ".join(set(sample_ids))))

        sample_types = [(sample.sample_type.decode('UTF-8') if isinstance(sample.sample_type, bytes) else sample.sample_type) for sample in samples]
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

    def join_fields(f, sep=", "):
        return sep.join([c.encode("utf-8") if isinstance(c, bytes) else c for c in f])


    seqr_fields['sex'] = join_fields(sorted(set(seqr_fields['sex'])))
    seqr_fields['genome (seqr)'] = join_fields(sorted(set(seqr_fields['genome (seqr)'])))
    seqr_fields['population (seqr)'] = join_fields(sorted(set(seqr_fields['population (seqr)'])))
    seqr_fields['sample id (seqr)'] = join_fields(sorted(set(seqr_fields['sample id (seqr)'])))
    seqr_fields['sample type (seqr)'] = join_fields(sorted(set(seqr_fields['sample type (seqr)'])))
    seqr_fields['analysis status (seqr)'] = join_fields(sorted(set(seqr_fields['analysis status (seqr)'])))

    seqr_fields['variant tags (seqr)'] = ", ".join(["%s %s tag(s)" % (n, key) for key, n in sorted(variant_tag_counters.items(), key=lambda x: x[0], reverse=True)]) + " "
    seqr_fields['variant notes (seqr)'] = "%s" % num_variant_notes

    seqr_fields['coded phenotype (seqr)'] = join_fields(seqr_fields['coded phenotype (seqr)']) + " "
    seqr_fields['anlaysis summary + notes (seqr)'] = join_fields(seqr_fields['anlaysis summary + notes (seqr)'], sep="\n") + " "
    seqr_fields['internal case review notes (seqr)'] = join_fields(seqr_fields['internal case review notes (seqr)'], sep="\n") + " "
    seqr_fields['cram path (seqr)'] = join_fields(sorted(set(seqr_fields['cram path (seqr)']))) + " "

seqr_fields_by_sample_id.values()

#sys.exit(0)
#%%
SEQR_INFO_COLUMNS = [
    'indiv (seqr)',
    'proj (seqr)',
    'fam (seqr)',
    'proj2 (seqr)',
    'fam2 (seqr)',
    'sex',
    'genome (seqr)',
    'population (seqr)',
    'sample id (seqr)',
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

final_df = final_df.merge(seqr_info_df, on="sample_id", how="left")

#%%
final_df

#%%

# data_paths_df left join Beryl's tables

beryls_supplementary_table_df = get_beryls_supplementary_table_df()
print("beryls_supplementary_table_df", beryls_supplementary_table_df.shape, beryls_supplementary_table_df.columns)

beryls_supplementary_table_df["Notes"] = beryls_supplementary_table_df["Notes on genetic diagnosis status"].fillna('')
for i, row in beryls_supplementary_table_df.iterrows():
    #print(i, ":", beryls_supplementary_table_df.at[i, "Notes"], ",", row["Notes from paper"])
    if row["Notes from paper"] and not isinstance(row["Notes from paper"], float):
        beryls_supplementary_table_df.at[i, "Notes"] += "\n" + row["Notes from paper"]

beryls_supplementary_table_df = beryls_supplementary_table_df[[
    "Sample ID",
    "Alias",
    "Clinical Diagnosis",
    "Sex",
    "Age at muscle biopsy",
    "Site of biopsy",
    "Notes",
    #"Previous NGS testing",
    #"Notes on genetic diagnosis status",
    #"Notes from paper",
]]

beryls_supplementary_table_df = beryls_supplementary_table_df.rename(columns={c: c.replace("\n", " ") + " (Beryl:Supp.)" for c in beryls_supplementary_table_df.columns if c != "Sample ID"})

final_df = left_join_beryls_table(final_df, beryls_supplementary_table_df)

#%%

beryls_rnaseq_probands_df = get_beryls_rnaseq_probands_df()
beryls_rnaseq_probands_df = beryls_rnaseq_probands_df.replace({
"Sample ID": {
    '62R_CaM (new)': '62R_CaM_3',
    '62R_CaM': '62R_CaM_2',
    'B13-07': 'B13-07-1M',
    '329FK_RR_R1': 'CLA_329FK_RR_2',
    '338FT_DM_R1': '338FT_DM_2',
    '214DF_AB_R1': 'CLA_214DF_AB_2',
    'B14-78.1M': 'B14-78-1-U',
}})

beryls_rnaseq_probands_df = beryls_rnaseq_probands_df[[
    "Sample ID",
    '%Contamin\n RNAseq',
    'Age at \nBiopsy',
    'Biopsy type',
    'Candidate \nVariants',
    'CanditateGenes\n(culprit,if solved)',
    'Collab PI',
    'Data details',
    'Data_type',
    'Ethnicity',
    'Include in manuscript?',
    "Phenotype", # "phenotype (Beryl)",
    'Short Phenotype',
    'Phenotype comments',
    'Other comments',
    'Variant \ntype',
    'Status',
    'Variant consequence',
]]

beryls_rnaseq_probands_df = beryls_rnaseq_probands_df.rename(columns={
    "Status": "Genetic diagnosis Status",
})

beryls_rnaseq_probands_df = beryls_rnaseq_probands_df.rename(columns={c: c.replace("\n", " ") + " (Beryl:Probands)" for c in beryls_rnaseq_probands_df.columns if c != "Sample ID"})

final_df = left_join_beryls_table(final_df, beryls_rnaseq_probands_df)


#%%

beryls_seqr_data_df = get_beryls_seqr_data_df()
beryls_seqr_data_df['Sample ID'] = beryls_seqr_data_df['Collaborator Participant ID']

beryls_seqr_data_df = beryls_seqr_data_df.replace({
"Sample ID": {
    'B09-24.1_UNKNOWN': 'B09-24-1RNA_UNKNOWN',
    'B13-07.1': 'B13-07-1M',
    '62R_CaM': '62R_CaM_2',
}})

for i, row in beryls_seqr_data_df.iterrows():
    #print(i, ":", beryls_supplementary_table_df.at[i, "Variant type(s)"], ",", row["Look at"])
    if row["Variant type(s)"] and not isinstance(row["Variant type(s)"], float):
        beryls_seqr_data_df.at[i, "Look at"] = (row["Variant type(s)"] or "") + ":" + (row["Look at"] or "")

beryls_seqr_data_df = beryls_seqr_data_df[[
    "Sample ID",
    'Phenotype',
    #'Variant type(s)',
    'Look at',
    'Notes',
    'RNA tissue (definitive source)',
]]

beryls_seqr_data_df = beryls_seqr_data_df.rename(columns={c: c.replace("\n", " ") + " (Beryl:Seqr-data)" for c in beryls_seqr_data_df.columns if c != "Sample ID"})

final_df = left_join_beryls_table(final_df, beryls_seqr_data_df)



#%%

df_export = final_df[[c for c in final_df.columns if c not in ('Sample ID', 'hg19_bam', 'rnaseqc_metrics')]]

df_export.iloc[270:280]

# %%


# export joined data to SEQR_INFO_AND_OTHER_METADATA_WORKSHEET
ws = get_seqr_info_and_other_metadata_worksheet()
set_with_dataframe(ws, df_export.fillna(''), resize=True)

print("Updated", ws.title)


# %%
