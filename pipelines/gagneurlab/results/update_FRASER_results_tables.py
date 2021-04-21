#%%

import glob
import os
import pandas as pd
from pprint import pprint
from collections import defaultdict
import hail as hl
hl.init(log="/dev/null")

from gagneurlab.gagneur_utils import get_rnaseq_truth_data_spreadsheet, get_RNASEQ_results_spreadsheet
from gspread_dataframe import set_with_dataframe
from gspread import WorksheetNotFound

#%%

# check results vs truth data

rows = get_rnaseq_truth_data_spreadsheet().worksheet("Events to detect").get()
truth_data_df = pd.DataFrame(data=rows[1:], columns=rows[0])

truth_data_df.columns

#%%

import re

truth_data_df[['chrom', 'start', 'end']] = truth_data_df["Junction (GRCh38)"].apply(lambda x: pd.Series(re.split("[:-]", str(x))))
truth_data_df['end'] = truth_data_df['end'].fillna(truth_data_df['start'].astype('int'))
truth_data_df.loc[:, 'chrom'] = truth_data_df.chrom.str.replace("chr", "")
truth_data_df.loc[:, 'start'] = truth_data_df.start.astype('int')
truth_data_df.loc[:, 'end'] = truth_data_df.end.astype('int')

#%%

star_combined_muscle_df = pd.read_table("~/project__rnaseq/code/rnaseq_methods/pipelines/gagneurlab//FRASER_results3/combined.muscle.153_samples.SJ.out.tsv.gz")
star_combined_muscle_df = star_combined_muscle_df[[
    "chrom", "start_1based", "end_1based", "strand", "intron_motif", "known_splice_junction", "unique_reads",
    "multi_mapped_reads", "maximum_overhang", "num_samples_with_this_junction", "num_samples_total",
]]

#star_combined_gtex_muscle_df = pd.read_table("~/project__rnaseq/code/rnaseq_methods/pipelines/gagneurlab/results/fraser/combined.gtex_muscle.803_samples.SJ.out.tsv.gz")
#star_combined_gtex_muscle_df = star_combined_gtex_muscle_df[["chrom", "start_1based", "end_1based", "strand", "intron_motif", "known_splice_junction", "unique_reads", "multi_mapped_reads", "maximum_overhang", "num_samples_with_this_junction", "num_samples_total"]]

#%%

for i, current_df in enumerate([star_combined_muscle_df, ]): # star_combined_gtex_muscle_df

    print(i, current_df.columns)
    current_df.loc[:, "chrom"] = current_df["chrom"].str.replace("chr", "")
    current_df.loc[:, "start_1based"] = current_df["start_1based"].astype('int')
    current_df.loc[:, "end_1based"] = current_df["end_1based"].astype('int')


#%%


#fraser_calls_df_with_GTEx = pd.read_table(
#    "~/project__rnaseq/code/rnaseq_methods/pipelines/gagneurlab/results/fraser/muscle/muscle_with_GTEX__253_samples_98B010A534__dpsi_0.1_padj_0.1_reads_2/"
#    "muscle_with_GTEX__253_samples_98B010A534_usingPCA_fds__psi5_q17__psi3_q14__psiSite_q20_padj_0.1__dpsi_0.1_results.tsv.gz").sort_values("padjust")

fraser_calls_df_without_GTEx = pd.read_table(
    "~/project__rnaseq/code/rnaseq_methods/pipelines/gagneurlab/FRASER_results3/muscle/muscle_without_GTEX__160_samples_D9B517A0F3/"
    "muscle_without_GTEX__160_samples_D9B517A0F3_usingPCA_fds__psi5_q5__psi3_q5__psiSite_q14_padj_0.1__dpsi_0.1_results.tsv.gz").sort_values("padjust")

#fraser_calls_df_without_GTEx_all_results = pd.read_table(
#    "~/project__rnaseq/code/rnaseq_methods/pipelines/gagneurlab/results/fraser/muscle_without_GTEX__153_samples_87843E873A/"
#    "muscle_without_GTEX__153_samples_87843E873A_usingPCA_fds__psi5_q5__psi3_q8__psiSite_q14_all_results.tsv.gz")
#fraser_calls_df_without_GTEx_all_results.set_index("Sample ID", inplace=True)

#outrider_calls_df_with_GTEx = pd.read_table(
#    "~/project__rnaseq/code/rnaseq_methods/pipelines/gagneurlab/results/outrider/muscle__153_samples_87843E873A_with_GTEX/"
#    "muscle__153_samples_87843E873A_with_GTEX__ods__q50_padj_0.05_results.tsv")

outrider_calls_df_without_GTEx = pd.read_table(
    "~/project__rnaseq/code/rnaseq_methods/pipelines/gagneurlab/OUTRIDER_results3/muscle__160_samples_D9B517A0F3_without_GTEX/"
    "muscle__160_samples_D9B517A0F3_without_GTEX__ods__q34_all_results.tsv"
)

#%%

all_rows_ranked_at_or_above_truth_junction = []
for with_GTEx, fraser_calls_df, outrider_calls_df in [
    (False, fraser_calls_df_without_GTEx, outrider_calls_df_without_GTEx)]: #,
    #(True, fraser_calls_df_with_GTEx, outrider_calls_df_with_GTEx)]:

    num_outrider_outliers_per_sample = []
    num_fraser_outliers_per_sample = []

    fraser_calls_df = fraser_calls_df.set_index(["seqnames", "start", "end"])
    fraser_calls_df["Is outlier in this # of samples"] = fraser_calls_df.groupby(["seqnames", "start", "end"])["sampleID"].nunique()
    fraser_calls_df = fraser_calls_df.reset_index()

    for _, truth_data_row in truth_data_df.iterrows():
        #if truth_data_row["Sample ID"] != "149BP_AB_M1":
        #    continue

        rows_ranked_at_or_above_truth_junction = []
        output_d = {
            "Truth Sample ID": truth_data_row["Sample ID"],
            "Truth Sample Alias": truth_data_row["Alias"],
            "Truth Junction": truth_data_row["Junction (GRCh38)"],
            "Truth Gene": truth_data_row["Gene"],
            "FRASER Included GTEx controls": with_GTEx,
        }

        truth_data_sample_id = str(truth_data_row["Sample ID"]).replace("MAAC031", "TOP_MAAC031_F_Muscle1").replace("49O_NM", "49O_NM_M1")

        print("---")
        print(f'Truth data: {truth_data_row["Sample ID"]}. {truth_data_row.chrom}:{truth_data_row.start}-{truth_data_row.end}')
        print_string = ""

        # check portcullis table:
        with hl.hadoop_open(f"gs://macarthurlab-rnaseq/batch_0/portcullis/{truth_data_sample_id}.portcullis_all.junctions.tab.gz", "r") as f:
            portcullis_all_df = pd.read_table(f)

        with hl.hadoop_open(f"gs://macarthurlab-rnaseq/batch_0/portcullis/{truth_data_sample_id}.portcullis_filtered.pass.junctions.tab.gz", "r") as f:
            portcullis_filtered_df = pd.read_table(f)

            # columns: ['index', 'refid', 'refname', 'reflen', 'start', 'end', 'size', 'left',
            #        'right', 'read-strand', 'ss-strand', 'consensus-strand', 'ss1', 'ss2',
            #        'canonical_ss', 'score', 'suspicious', 'pfp', 'nb_raw_aln',
            #        'nb_dist_aln', 'nb_us_aln', 'nb_ms_aln', 'nb_um_aln', 'nb_mm_aln',
            #        'nb_bpp_aln', 'nb_ppp_aln', 'nb_rel_aln', 'rel2raw', 'nb_r1_pos',
            #        'nb_r1_neg', 'nb_r2_pos', 'nb_r2_neg', 'entropy', 'mean_mismatches',
            #        'mean_readlen', 'max_min_anc', 'maxmmes', 'intron_score', 'hamming5p',
            #        'hamming3p', 'coding', 'pws', 'splice_sig', 'uniq_junc', 'primary_junc',
            #        'nb_up_juncs', 'nb_down_juncs', 'dist_2_up_junc', 'dist_2_down_junc',
            #        'dist_nearest_junc', 'mm_score', 'coverage', 'up_aln', 'down_aln',
            #        'nb_samples', 'JAD01', 'JAD02', 'JAD03', 'JAD04', 'JAD05', 'JAD06',
            #        'JAD07', 'JAD08', 'JAD09', 'JAD10', 'JAD11', 'JAD12', 'JAD13', 'JAD14',
            #        'JAD15', 'JAD16', 'JAD17', 'JAD18', 'JAD19', 'JAD20']
        portcullis_tables = {}
        for portcullis_df_label, portcullis_df in [('portcullis_passed_filter', portcullis_filtered_df), ("portcullis_all", portcullis_all_df)]:
            portcullis_df.loc[:, 'chrom'] = portcullis_df.refname.str.replace("chr", "")
            portcullis_df.loc[:, 'start'] = portcullis_df.start.astype('int')
            portcullis_df.loc[:, 'end'] = portcullis_df.end.astype('int')
            portcullis_df = portcullis_df.set_index(["chrom", "start", "end"], drop=False)
            portcullis_tables[portcullis_df_label] = portcullis_df

        for portcullis_df_label, portcullis_df in portcullis_tables.items():
            try:
                portcullis_df_rows = [portcullis_df.loc[(truth_data_row.chrom, truth_data_row.start, truth_data_row.end - 1), :]]
            except Exception as e:
                print(e)
                portcullis_df_rows = []

            #portcullis_df_rows = portcullis_df[(
            #    (portcullis_df.chrom == truth_data_row.chrom) &
            #    (portcullis_df.start == truth_data_row.start) &
            #    (portcullis_df.end == truth_data_row.end - 1))]
            print(portcullis_df_label, len(portcullis_df_rows))
            if len(portcullis_df_rows) > 0:
                output_d[f"portcullis_table"] = portcullis_df_label
                for key in "score", "suspicious", "mean_mismatches", "entropy", "max_min_anc", "maxmmes", "splice_sig", "pws":
                    output_d[f"portcullis_{key}"] = portcullis_df_rows[0][key]
                print(output_d)
                break

        sample_fraser_calls_df = fraser_calls_df[fraser_calls_df["sampleID"] == truth_data_row["Sample ID"]]
        num_fraser_outliers_per_sample.append(len(sample_fraser_calls_df))

        sample_fraser_calls_df.loc[:, "chrom"] = sample_fraser_calls_df.seqnames.str.replace("chr", "")
        sample_fraser_calls_df.loc[:, 'start'] = sample_fraser_calls_df.start.astype('int')
        sample_fraser_calls_df.loc[:, 'end'] = sample_fraser_calls_df.end.astype('int')

        truth_idx = None
        # if truth_idx is None:
        #     sample_fraser_junction_match = sample_fraser_calls_df[
        #         (sample_fraser_calls_df.chrom == truth_data_row.chrom) &
        #         (sample_fraser_calls_df.start - 1 == truth_data_row.start) &
        #         (sample_fraser_calls_df.end == truth_data_row.end)]
        #     if len(sample_fraser_junction_match) > 0:
        #         print("### Exact match junction")
        #         truth_idx = sample_fraser_junction_match.index[0]

        #if truth_idx is None:
        #    sample_fraser_junction_match = sample_fraser_calls_df[
        #        (sample_fraser_calls_df.chrom == truth_data_row.chrom) & (
        #                (abs(sample_fraser_calls_df.start - truth_data_row.start) <= 1) & (abs(sample_fraser_calls_df.end - truth_data_row.end) <= 1)
        #        )]
        #    if len(sample_fraser_junction_match) > 0:
        #        print("### Exact match between truth junction & fraser results table")
        #        truth_idx = sample_fraser_junction_match.index[0]
        #        output_d["FRASER junction match type"] = "both junction coords"

        # if truth_idx is None:
        #     sample_fraser_junction_match = sample_fraser_calls_df[
        #         ((sample_fraser_calls_df.chrom == truth_data_row.chrom) & (sample_fraser_calls_df.start == truth_data_row.start)) |
        #         ((sample_fraser_calls_df.chrom == truth_data_row.chrom) & (sample_fraser_calls_df.end == truth_data_row.end))]
        #     if len(sample_fraser_junction_match) > 0:
        #         print("### Exact match position")
        #         truth_idx = sample_fraser_junction_match.index[0]

        if truth_idx is None:
            sample_fraser_junction_match = sample_fraser_calls_df[
                (sample_fraser_calls_df.chrom == truth_data_row.chrom) & (
                        (abs(sample_fraser_calls_df.start - truth_data_row.start) <= 1) | (abs(sample_fraser_calls_df.end - truth_data_row.end) <= 1)
                )]
            if len(sample_fraser_junction_match) > 0:
                print("#### Inexact match position between truth junction & fraser results table")
                truth_idx = sample_fraser_junction_match.index[0]
                output_d["FRASER junction match type"] = "one junction coord"

        if truth_idx is None:
            sample_fraser_junction_match = sample_fraser_calls_df[sample_fraser_calls_df.hgncSymbol == truth_data_row.Gene]
            if len(sample_fraser_junction_match) > 0:
                print("#### Gene match position between truth junction & fraser results table")
                truth_idx = sample_fraser_junction_match.index[0]
                output_d["FRASER junction match type"] = "gene"

        if truth_idx is None:
            list_of_fraser_table_rows = []
            print("#### No match found in fraser results table")
            output_d["FRASER junction match type"] = "no match"
        else:
            list_of_fraser_table_rows = list(sample_fraser_calls_df.loc[:truth_idx,].iterrows())

        print(f"## FRASER Truth Rank: {len(list_of_fraser_table_rows)}")

        output_d["FRASER num junctions ranked at or above truth junction"] = len(list_of_fraser_table_rows)
        output_d["FRASER num junctions significant in this sample"] = len(sample_fraser_calls_df)

        if len(list_of_fraser_table_rows) == 0:
            list_of_fraser_table_rows = [
                (truth_idx, {"chrom": truth_data_row.chrom, "start": truth_data_row.start, "end": truth_data_row.end})
            ]

        # check  OUTRIDER
        sample_outrider_calls_df = outrider_calls_df[outrider_calls_df["sampleID"] == truth_data_row["Sample ID"]]
        output_d["OUTRIDER num junctions significant in this sample"] = len(sample_outrider_calls_df)
        num_outrider_outliers_per_sample.append(output_d["OUTRIDER num junctions significant in this sample"])

        sample_outrider_gene_match = sample_outrider_calls_df[sample_outrider_calls_df.geneID == truth_data_row.Gene]
        outrider_truth_idx = None
        if outrider_truth_idx is None:
            list_of_outrider_table_rows = []
            print("#### No match found in outrider results table")
        else:
            list_of_outrider_table_rows = list(sample_outrider_gene_match.loc[:outrider_truth_idx,].iterrows())

        print(f"OUTRIDER truth rank? {len(list_of_outrider_table_rows)} out of {len(sample_outrider_calls_df)}")


        with hl.hadoop_open(f'gs://macarthurlab-rnaseq/batch_0/star/{truth_data_sample_id}.SJ.out.tab.gz', "r") as f:
            star_df = pd.read_table(f, names=["chrom", "start_1based", "end_1based", "strand", "intron_motif", "known_splice_junction", "unique_reads", "multi_mapped_reads", "maximum_overhang"])
            star_df.loc[:,"chrom"] = star_df["chrom"].str.replace("chr", "")
            star_df.loc[:,"start_1based"] = star_df["start_1based"].astype('int')
            star_df.loc[:,"end_1based"] = star_df["end_1based"].astype('int')

        for i, (idx, row) in enumerate(list_of_fraser_table_rows):
            output_d_copy = dict(output_d)
            output_d_copy['is_truth_junction'] = idx == truth_idx
            #if not output_d_copy['is_truth_junction']:
            #    continue
            if not isinstance(row,  dict):
                output_d_copy.update(row.to_dict())
                output_d_copy[f"FRASER Junction Gene"] = row["hgncSymbol"]
                output_d_copy["FRASER Junction"] = f"{row.seqnames}:{row.start}-{row.end}"
                output_d_copy[f"FRASER Metric Type"] = row["type"]
                output_d_copy[f"FRASER deltaPsi"] = row["deltaPsi"]
                output_d_copy[f"FRASER Padj"] = row["padjust"]
                output_d_copy[f"FRASER zScore"] = row["zScore"]
                output_d_copy[f"FRASER Junction counts"] = row["counts"]
                output_d_copy[f"FRASER Junction total counts"] = row["totalCounts"]
                output_d_copy['rank'] = i + 1
            else:
                print("Processing truth data with no matching FRASER records")
                output_d_copy["FRASER Junction Gene"] = ""
                output_d_copy["FRASER Junction"] = ""
                output_d_copy['rank'] = 0

            #all_rows_ranked_at_or_above_truth_junction.append(output_d_copy)
            #continue

            # check portcullis tables
            for portcullis_df_label, portcullis_df in portcullis_tables.items():

                portcullis_df_rows = [r for _, r in portcullis_df[(
                    (portcullis_df.chrom == row['chrom']) &
                    (abs(portcullis_df.start - int(row['start'])) <= 1) &
                    (abs(portcullis_df.end - int(row['end'])) <= 1)
                )].iterrows()]

                print(portcullis_df_label, len(portcullis_df_rows))
                if len(portcullis_df_rows) > 0:
                    output_d_copy[f"portcullis_table"] = portcullis_df_label
                    for key in "score", "suspicious", "mean_mismatches", "entropy", "max_min_anc", "maxmmes", "splice_sig", "pws":
                        output_d_copy[f"portcullis_{key}"] = portcullis_df_rows[0][key]
                    #print(output_d_copy)
                    break

            matching_df_counter = 0
            for key_label, current_df in [
                ("STAR_sample", star_df),
                ("STAR_combined_muscle", star_combined_muscle_df),
                #("STAR_combined_GTEx_muscle", star_combined_gtex_muscle_df),
            ]:
                current_df_matching_rows = current_df[(
                    (current_df["chrom"].replace("chr", "") == row['chrom'].replace("chr", "")) &
                    (abs(current_df["start_1based"] - int(row['start'])) <= 1) &
                    (abs(current_df["end_1based"] - int(row['end'])) <= 1)
                )]

                for key in "unique_reads", "multi_mapped_reads":
                    output_d_copy[f"{key_label}_{key}"] = 0

                if len(current_df_matching_rows) == 1:
                    output_d_copy[f'FRASER Junction Has 0 matches in {key_label}'] = False
                    matching_df_counter += len(current_df_matching_rows)
                    print(f'    {truth_data_row["Sample ID"]} ({i}): found exact match between fraser row {row["chrom"]}:{row["start"]}-{row["end"]} and {key_label} row {current_df_matching_rows["chrom"].iloc[0]}:{current_df_matching_rows["start_1based"].iloc[0]}-{current_df_matching_rows["end_1based"].iloc[0]}')
                    for key, value in current_df_matching_rows.iloc[0].to_dict().items():
                        if key not in {"chrom", "start_1based", "end_1based", "strand", "intron_motif"}:
                            output_d_copy[f"{key_label}_{key}"] = value
                            #if truth_data_sample_id == "149BP_AB_M1":
                            #    print(f"{key_label}_{key} : {value}")
                elif len(current_df_matching_rows) > 1:
                    output_d_copy[f'FRASER Junction Has 0 matches in {key_label}'] = False
                    print(f'    ==== {truth_data_row["Sample ID"]} ({i}): fraser row {row.chrom}:{row.start}-{row.end} matches too many rows ({len(current_df_matching_rows)} rows) in {key_label}. ')
                else:
                    output_d_copy[f'FRASER Junction Has 0 matches in {key_label}'] = True
                    if not isinstance(row,  dict) and (output_d_copy[f"FRASER Metric Type"] == "psi5" or output_d_copy[f"FRASER Metric Type"] == "psi3"):
                        print(f'    #### {truth_data_row["Sample ID"]} ({i}): fraser row {row["chrom"]}:{row["start"]}-{row["end"]}  {output_d_copy[f"FRASER Metric Type"]} matches zero rows in {key_label}. ')
                    else:
                        print(f'    {truth_data_row["Sample ID"]} ({i}): fraser row {row["chrom"]}:{row["start"]}-{row["end"]}  matches matches zero rows in {key_label}')


            #print(output_d_copy['seqnames'], output_d_copy['start'], output_d_copy['end'], " matching star rows ", matching_df_counter)
            rows_ranked_at_or_above_truth_junction.append(output_d_copy)

        for d in rows_ranked_at_or_above_truth_junction:
            d["normalized rank"] = d["rank"] / len(rows_ranked_at_or_above_truth_junction)

        all_rows_ranked_at_or_above_truth_junction.extend(rows_ranked_at_or_above_truth_junction)

        print(f'{truth_data_row["Sample ID"]} ({truth_data_row["Alias"]}): {i+1} results processed. Event {output_d_copy["FRASER Junction"]}. Truth Junction: {truth_data_row["Junction (GRCh38)"]} \n')
        print("")
        #if output_d["Detected"] != "Not Detected":
        #    print(print_string)
    print("FRASER calls per sample", list(sorted(num_fraser_outliers_per_sample)))
    print("OUTRIDER calls per sample", list(sorted(num_outrider_outliers_per_sample)))

pd.DataFrame(all_rows_ranked_at_or_above_truth_junction).to_csv("all_rows_ranked_at_or_above_truth_junction.tsv", header=True, index=False, sep="\t")



#%%

all_fraser_calls_df_without_GTEx = pd.read_table(
    "~/project__rnaseq/code/rnaseq_methods/pipelines/gagneurlab/results/fraser/muscle_without_GTEX__153_samples_87843E873A/"
    "muscle_without_GTEX__153_samples_87843E873A_usingPCA_fds__psi5_q5__psi3_q8__psiSite_q14_all_results.tsv.gz")


#%%

all_fraser_calls_df_without_GTEx.columns


#%%
#%%
#%%
# gsutil -m cp gs://macarthurlab-rnaseq/gagneur/outrider/results/*with*.tsv.gz .

all_tables = defaultdict(list)
results_tables = "/Users/weisburd/project__rnaseq/code/rnaseq_methods/pipelines/gagneurlab/FRASER_results1/using_AE/*.tsv.gz"
for path in [p for p in glob.glob(results_tables)]:
    name = os.path.basename(path)
    label = name.split("_with")[0]
    if "_without_GTEX" in name:
        has_GTEX = False
    elif "_with_GTEX" in name:
        has_GTEX = True
    else:
        raise ValueError(f"Unexpected filename: {name}")
    all_tables[label].append({'has_GTEX': has_GTEX, 'path': path})

pprint(sorted(all_tables))

#%%


def read_table(path, padj_threshold=0.10):
    df = pd.read_table(path)
    df = df[df.padjust < padj_threshold]
    df['geneID'] = df['geneID'].apply(lambda s: s.split(".")[0]) # remove gene versions
    df = df.set_index(["sampleID", "geneID"])
    return df


results = {}
for label, tables in all_tables.items():
    if len(tables) > 2:
        raise ValueError(f"More than 2 tables found for: {label}")

    if len(tables) == 1:
        results[label] = read_table(tables[0]['path']).reset_index()
        continue
    tables.sort(key=lambda x: x['has_GTEX'])
    t1 = read_table(tables[0]['path'])
    t2 = read_table(tables[1]['path'])

    #print(t1.columns)  # 'sampleID', 'geneID', 'pValue', 'padjust', 'zScore', 'rawcounts', 'q'
    result = t1.join(t2,
        how="outer",
        lsuffix=("_with_GTEX" if tables[0]['has_GTEX'] else "_without_GTEX"),
        rsuffix=("_with_GTEX" if tables[1]['has_GTEX'] else "_without_GTEX")).reset_index()
    results[label] = result


#%%

# https://docs.google.com/spreadsheets/d/1mvUcRDyAINf1utpF4TqplE4lAS_8ibglZA7iDFx7k2o/edit?usp=sharing

spreadsheet = get_rnaseq_truth_data_spreadsheet()

for label, results_df in results.items():
    worksheet_name = f"{label} (auto)"
    try:
        worksheet = spreadsheet.worksheet(worksheet_name)
    except WorksheetNotFound:
        worksheet = spreadsheet.add_worksheet(worksheet_name, 1, 1)
        print("Created", worksheet.title)

    set_with_dataframe(worksheet, results_df.reset_index().fillna(''), resize=True)
    print("Updated", worksheet.title)

#%%

# use OMIM API to add OMIM info


#%%

# seqr gene lists


#%%

# RNA-seq sample metadata

#%%

# add constraint info
constraint_df = pd.read_table("https://storage.googleapis.com/gnomad-public/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz", compression="gzip")

#for c in sorted(constraint_df.columns):
#    print(c)

for k in results.keys():
    result = results[k]
    result.geneID = result.geneID.apply(lambda s: s.split(".")[0])
    break

#%%

# add gene list column



#%%
# Uplaod tables to Google Sheets

spreadsheet = get_RNASEQ_results_spreadsheet()

#%%

# upload tables
results_tsvs = {
    #"whole_blood": "~/project__rnaseq/code/rnaseq_methods/pipelines/gagneurlab/results/fraser/fibroblasts/fibroblasts_with_GTEX__186_samples_5A9664E950__dpsi_0.1_padj_0.1/fibroblasts_with_GTEX__186_samples_5A9664E950_usingPCA_fds__psi5_q8__psi3_q8__psiSite_q11_padj_0.1__deltapsi_0.1_results.tsv",

    #"FRASER: muscle without GTEx (153 samples, padj_0.1, dpsi_0.1)": "~/project__rnaseq/code/rnaseq_methods/pipelines/gagneurlab/results/fraser/muscle/muscle_without_GTEX__153_samples_87843E873A/muscle_without_GTEX__153_samples_87843E873A_usingPCA_fds__psi5_q5__psi3_q8__psiSite_q14_padj_0.1__dpsi_0.1_results.tsv",
    #"FRASER: fibroblasts without GTEx (86 samples, padj_0.1, dpsi_0.1)": "~/project__rnaseq/code/rnaseq_methods/pipelines/gagneurlab/results/fraser/fibroblasts/fibroblasts_without_GTEX__86_samples_867E0896FE/fibroblasts_without_GTEX__86_samples_867E0896FE_usingPCA_fds__psi5_q2__psi3_q2__psiSite_q5_padj_0.1__dpsi_0.1_results.tsv",
    #"FRASER: whole blood with GTEx (107 samples, padj_0.1, dpsi_0.1)": "~/project__rnaseq/code/rnaseq_methods/pipelines/gagneurlab/results/fraser/whole_blood/whole_blood_with_GTEX__107_samples_F3E00680DC/whole_blood_with_GTEX__107_samples_F3E00680DC_usingPCA_fds__psi5_q8__psi3_q8__psiSite_q8_padj_0.1__dpsi_0.1_results.tsv.gz",

    #"OUTRIDER: muscle without GTEx (153 samples, padj0.05)": "~/project__rnaseq/code/rnaseq_methods/pipelines/gagneurlab/results/outrider/muscle__153_samples_87843E873A_without_GTEX/muscle__153_samples_87843E873A_without_GTEX__ods__q32_padj_0.05_results.tsv",
    #"OUTRIDER: fibroblasts without GTEx (86 samples, padj0.05)": "~/project__rnaseq/code/rnaseq_methods/pipelines/gagneurlab/results/outrider/fibroblasts__86_samples_867E0896FE_without_GTEX/fibroblasts__86_samples_867E0896FE_without_GTEX__ods__q18_padj_0.05_results.tsv",
    #"OUTRIDER: whole blood with GTEx (107 samples, padj0.05)": "~/project__rnaseq/code/rnaseq_methods/pipelines/gagneurlab/results/outrider/whole_blood__7_samples_0BFE4A559F_with_GTEX/whole_blood__7_samples_0BFE4A559F_with_GTEX__ods__q22_padj_0.05_results.tsv",
}


for label, results_path in results_tsvs.items():
    worksheet_name = f"{label} (auto)"
    try:
        worksheet = spreadsheet.worksheet(worksheet_name)
    except WorksheetNotFound:
        worksheet = spreadsheet.add_worksheet(worksheet_name, 1, 1)
        print("Created", worksheet.title)

    results_df = pd.read_table(results_path)
    set_with_dataframe(worksheet, results_df.reset_index().fillna(''), resize=True)
    print("Updated", worksheet.title)
