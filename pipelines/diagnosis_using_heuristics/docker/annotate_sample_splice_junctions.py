
"""
TODO

take individual junctions

left-join with:
    portcullis
    all 4 GTEx tissues  (# of samples with junction, # of samples total, min, max, median read support per sample)
    variants
    OMIM, constraint, gene lists, etc.

"""
# TODO OUTRIDER data
# TODO variants from sample

import glob
import os
import pandas as pd

print(os.getcwd())
os.chdir("/Users/weisburd/project__rnaseq/code/rnaseq_methods/pipelines/diagnosis_using_heuristics/docker/test")

#%%

header_fields = "sample_id	alias	gene	junction".split("\t")

# B11-48-1M	C1	POMGNT1	chr1:46189457
# 253DY_HA_M1	C4	DMD	chrX:32518008

rows = [dict(zip(header_fields, row.split("\t"))) for row in """
126BG_CB_M1	D1	LARGE1	22:33337801-33919994
41M_MW_M1	D2	PYROXD1	12:21440448-21449562
153BR_JB_M1	D7	NEB	2:151672680-151675286
163BV_JE_M1	D8	NEB	2:151680828-151684777
210DB_BW_M1	D9	RYR1	19:38440864-38443557
T1244	D10	TTN	2:178622008-178624464
361AL_CC_M1	D11	COL6A3	2:237358583-237359205
149BP_AB_M1	N25	NEB	2:151498552-151499297
UC84-1RNA	E1	NEB	2:151531896-151534214
247DT_SH_M1	E2	NEB	2:151687733-151688371 
373HQ_BTG_M1	E4	TTN	2:178580609-178581905
251DW_SD_M1	C3	DMD	X:31729748-31819974
MAAC031	C7	DMD	X:31590401-31595570
26I_SK_M1	C9	NEB	2:151531896-151533382
CLA_180CJ_DP_2	C11	RYR1	19:38467812-38468965
UC316-1M	N31	COL6A1	21:45989778-45989893
UC393-1M	N32	COL6A1	21:45989778-45989893
49O_NM	N33	DMD	X:32256704-32287528
""".strip().split("\n")]


truth_junctions = pd.DataFrame(rows, columns=header_fields)

truth_junctions

def parse_interval(i):

    try:
        chrom, start_end = i.split(":")
        start, end = map(int, start_end.replace(",", "").split("-"))
        return chrom, start, end
    except Exception as e:
        raise ValueError(f"Couldn't parse interval: {i}. {e}")


def add_chrom_start_end_columns(row):
    try:
        row["chrom"], row["start_1based"], row["end_1based"] = parse_interval(row["junction"])
        row["chrom"] = "chr" + row["chrom"].replace("chr", "")  # convert 0-based to 1-based
        row["start_1based"] += 1  # convert 0-based to 1-based
    except Exception as e:
        row["chrom"] = row["start_1based"] = row["end_1based"] = None
        print(f"Skipping {row.sample_id}: {e}")

    print(row.to_dict())
    return row


truth_junctions = truth_junctions.apply(add_chrom_start_end_columns, axis=1)

truth_junctions.set_index("sample_id", inplace=True)

#truth_junctions

for path in glob.glob("./vcf/*.vcf.bgz"):
    sample_name = os.path.basename(path).split(".")[0]
    if sample_name in set(truth_junctions.index):
        truth_junctions.loc[sample_name, "vcf_path"] = path


for path in glob.glob("./junctions/*.SJ.out.tab.gz"):
    sample_name = os.path.basename(path).split(".")[0]
    sample_name = sample_name.replace("TOP_MAAC031_F_Muscle1", "MAAC031")
    sample_name = sample_name.replace("49O_NM_M1", "49O_NM")

    if sample_name in set(truth_junctions.index):
        truth_junctions.loc[sample_name, "junctions_path"] = path
    else:
        print(sample_name)


for path in glob.glob("./portcullis_all/*.tab"):
    sample_name = os.path.basename(path).split(".")[0]
    sample_name = sample_name.replace("TOP_MAAC031_F_Muscle1", "MAAC031")
    sample_name = sample_name.replace("49O_NM_M1", "49O_NM")

    if sample_name in set(truth_junctions.index):
        truth_junctions.loc[sample_name, "portcullis_all_path"] = path
    else:
        print(sample_name)


for path in glob.glob("./portcullis_filtered/*.tab"):
    sample_name = os.path.basename(path).split(".")[0]
    sample_name = sample_name.replace("TOP_MAAC031_F_Muscle1", "MAAC031")
    sample_name = sample_name.replace("49O_NM_M1", "49O_NM")

    if sample_name in set(truth_junctions.index):
        truth_junctions.loc[sample_name, "portcullis_filtered_path"] = path
    else:
        print(sample_name)

truth_junctions = truth_junctions.reset_index()

truth_junctions.set_index(["chrom", "start_1based", "end_1based"], inplace=True)

truth_junctions


#%%
combined_muscle_samples_df = pd.read_table("./combined_tables/combined.muscle.162_samples.SJ.out.tsv.gz")
combined_muscle_gtex_df = pd.read_table("./combined_tables/gtex/combined.muscle.100_gtex_only_samples.SJ.out.tsv.gz")

combined_muscle_samples_df.set_index(["chrom", "start_1based", "end_1based"], inplace=True)
combined_muscle_gtex_df.set_index(["chrom", "start_1based", "end_1based"], inplace=True)

#%%
combined_muscle_samples_df.columns = combined_muscle_samples_df.columns.map(lambda x: str(x) + '_combined_samples')
combined_muscle_gtex_df.columns = combined_muscle_gtex_df.columns.map(lambda x: str(x) + '_combined_gtex')

#%%

truth_junctions = truth_junctions.join(combined_muscle_samples_df, how="left")
truth_junctions = truth_junctions.join(combined_muscle_gtex_df, how="left")

truth_junctions.drop(columns=[
    f"{c}{s}" for c in ["strand", "intron_motif", "strand_counter", "known_splice_junction"] for s in ["_combined_gtex", "_combined_samples"]
], inplace=True)

#%%


for truth_junctions