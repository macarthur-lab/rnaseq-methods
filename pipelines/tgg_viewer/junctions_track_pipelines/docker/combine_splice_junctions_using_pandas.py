import argparse
from collections import defaultdict
import gc
import itertools
import logging
import pandas as pd
import numpy as np
import psutil
import os
import re

#pd.options.mode.chained_assignment = 'raise'

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s: %(message)s')

COLUMN_TYPES = {
    'chrom': 'str',
    'start_1based': 'int32',
    'end_1based': 'int32',
    'strand': 'int8',
    'intron_motif': 'int8',
    'known_splice_junction': 'int8',
    'unique_reads': 'int32',
    'multi_mapped_reads': 'int32',
    'total_reads': 'int32',
    'maximum_overhang': 'int16',
    'num_samples_with_this_junction': 'int32',
    'num_samples_total': 'int32',
    'max_per_sample_unique_reads': 'int32',
    'max_per_sample_total_reads': 'int32',
    'strand_counter': 'int32',
}

COLUMN_NAMES = [
    'strand',
    'intron_motif',
    'known_splice_junction',
    'unique_reads',
    'multi_mapped_reads',
    'maximum_overhang',
    'total_reads',
    'num_samples_with_this_junction',
    'num_samples_total',
    'max_per_sample_unique_reads',
    'max_per_sample_total_reads',
    'strand_counter',
    'sample_id',
]

AGGREGATED_COLUMNS_TO_ADD_TO_INDIVIDUAL_TABLES = [
    'num_samples_with_this_junction',
    'num_samples_total',
    'max_per_sample_unique_reads',
    'max_per_sample_total_reads',
    'sample_id',
]

def parse_args():
    p = argparse.ArgumentParser(description="Combines multiple SJ.out.tab tables into one combined .bed file which can be loaded into TGG-viewer to visualize splice junctions")
    p.add_argument("-n", "--batch-size", type=int, default=50, help="How many tables to merge in memory before writing to disk. Bigger batch sizes use more memory.")
    p.add_argument("-o", "--output-path", help="Combined .bed output path")
    p.add_argument("--normalize-read-counts", action="store_true", help="whether to normalize unique- and multi-mapped read counts rather than just summing them across input tables")
    p.add_argument("--save-individual-tables", action="store_true", help="Also export individual .bed files with additional columns")
    p.add_argument("--discard-sample-id-column", action="store_true", help="Don't add a sample_id column with a list of samples that have each splice junction. This increases compute time by 10x.")
    p.add_argument("paths", nargs="+", help="Paths of 1 or more input SJ.out.tab tables")
    args = p.parse_args()

    if args.normalize_read_counts:
        for column_name in 'unique_reads', 'multi_mapped_reads', 'total_reads', 'max_per_sample_unique_reads', 'max_per_sample_total_reads':
            COLUMN_TYPES[column_name] = 'float32'

    return args


def batched_iter(iterable, batch_size=1):
    it = iter(iterable)
    while True:
        batch = tuple(itertools.islice(it, batch_size))
        if not batch:
            break
        yield batch


prev_memory_bytes = 0
def print_memory_stats(message="", run_gc=False):
    global prev_memory_bytes
    if message:
        message = " - " + message
    if run_gc:
        gc.collect()

    memory_bytes = psutil.Process(os.getpid()).memory_info().rss

    logging.info(f'memory used {message}: {memory_bytes//10**6} Mb    delta: {(memory_bytes - prev_memory_bytes)//10**6} Mb')
    prev_memory_bytes = memory_bytes


def read_SJ_out_tab(path, i=None):
    df = pd.read_csv(
        path,
        names=[
            'chrom', 'start_1based', 'end_1based',
            f'strand',
            f'intron_motif',
            f'known_splice_junction',
            f'unique_reads',
            f'multi_mapped_reads',
            f'maximum_overhang'],
        index_col=['chrom', 'start_1based', 'end_1based'],
        dtype=COLUMN_TYPES,
        sep='\t')

    df.loc[:, 'total_reads'] = df["unique_reads"] + df["multi_mapped_reads"]

    df.loc[:, 'num_samples_with_this_junction'] = np.int32(1)
    df.loc[:, 'num_samples_total'] = np.int32(1)
    df.loc[:, 'max_per_sample_unique_reads'] = np.int32(0)
    df.loc[:, 'max_per_sample_total_reads'] = np.int32(0)

    df.loc[:, 'strand_counter'] = df['strand'].apply(lambda s: 1 if s == 1 else (-1 if s == 2 else 0)).astype('int32')
    df.loc[:, 'sample_id'] = os.path.basename(path).replace(".SJ.out.tab", "").replace(".gz", "")

    # print some stats
    print(f"   {df.unique_reads.sum()/1_000_000:0.1f} million uniquely-mapped reads")
    print(f"   {100*df.multi_mapped_reads.sum()/float(df.unique_reads.sum() + df.multi_mapped_reads.sum()):0.0f}% multi-mapped")

    if i is not None:
        name_mapping = dict(zip(df.columns, [f"{c}_{i}" for c in df.columns]))
        df.rename(columns=name_mapping, inplace=True)

    return df


"""
How the normalization works:

2 samples: 
800 reads,  1200 reads
total_unique_reads_across_all_samples = 1000
scalars:
1.25 (5/4), 0.83  (5/6)

junction1: 10,  15   =>  10*1.25/2  + 15*0.83/2 =>  normalized count: 12.475   
junction2: 10,  5    =>  10*1.25/2  +  5*0.83/2 =>  normalized count:  8.325 
"""


def process_all_batches(args):
    total_samples = len(args.paths)
    logging.info(f'Processing {total_samples} tables')

    if args.normalize_read_counts:
        # read all the tables to compute total_unique_reads_across_all_samples and average_unique_reads_per_sample
        total_unique_reads_across_all_samples = 0
        for path in args.paths:
            df = read_SJ_out_tab(path)
            if args.discard_sample_id_column:
                df.drop(columns=['sample_id'], inplace=True)
            total_unique_reads_across_all_samples += df.unique_reads.sum()

        average_unique_reads_per_sample = total_unique_reads_across_all_samples/float(total_samples)

    all_combined_SJ_out_tab_df = pd.DataFrame()
    sample_i = 0
    for batch_number, batch_sample_paths in enumerate(batched_iter(args.paths, args.batch_size)):
        tables_in_batch = []
        batch_start_i = sample_i
        for path in batch_sample_paths:
            logging.info(f"Batch {batch_number}, Table {sample_i}: {path}")
            current_SJ_out_tab_df = read_SJ_out_tab(path, sample_i)

            if args.normalize_read_counts:
                unique_reads_in_sample = current_SJ_out_tab_df[f"unique_reads_{sample_i}"].sum()
                scalar = average_unique_reads_per_sample / float(unique_reads_in_sample)

                current_SJ_out_tab_df.loc[:, f"unique_reads_{sample_i}"] *= scalar / float(total_samples)
                current_SJ_out_tab_df.loc[:, f"multi_mapped_reads_{sample_i}"] *= scalar / float(total_samples)
                current_SJ_out_tab_df.loc[:, f"total_reads_{sample_i}"] *= scalar / float(total_samples)

                current_SJ_out_tab_df.loc[:, f"unique_reads_{sample_i}"] = current_SJ_out_tab_df[f"unique_reads_{sample_i}"].round(decimals=3)
                current_SJ_out_tab_df.loc[:, f"multi_mapped_reads_{sample_i}"] = current_SJ_out_tab_df[f"multi_mapped_reads_{sample_i}"].round(decimals=3)
                current_SJ_out_tab_df.loc[:, f"total_reads_{sample_i}"] = current_SJ_out_tab_df[f"total_reads_{sample_i}"].round(decimals=3)
                logging.info(f"{path} has {int(unique_reads_in_sample)} total unique reads while the "
                             f"per-sample average is {average_unique_reads_per_sample}. Scaling read counts "
                             f"by {scalar} and dividing by {total_samples}")


            tables_in_batch.append(current_SJ_out_tab_df)
            sample_i += 1
        batch_end_i = sample_i

        all_combined_SJ_out_tab_df = all_combined_SJ_out_tab_df.join(tables_in_batch, how="outer")

        current_batch_columns = {}
        for column in COLUMN_NAMES:
            current_batch_columns[column] = [f'{column}_{k}' for k in range(batch_start_i, batch_end_i)]
            if batch_number > 0: current_batch_columns[column].append(column)

        # merge intron_motif columns across samples in the current batch
        intron_motif_columns_df = all_combined_SJ_out_tab_df[current_batch_columns['intron_motif']]
        intron_motif_columns_df.ffill(axis=1, inplace=True)
        intron_motif_columns_df = intron_motif_columns_df.iloc[:,-1].astype(COLUMN_TYPES['intron_motif'])
        all_combined_SJ_out_tab_df.loc[:, 'intron_motif'] = intron_motif_columns_df

        # merge known_splice_junction columns across samples in the current batch
        known_splice_junction_df = all_combined_SJ_out_tab_df[current_batch_columns['known_splice_junction']]
        known_splice_junction_df.replace(0, np.nan, inplace=True)
        known_splice_junction_df.ffill(axis=1, inplace=True)
        known_splice_junction_df = known_splice_junction_df.iloc[:,-1]
        known_splice_junction_df.fillna(0, inplace=True)
        known_splice_junction_df = known_splice_junction_df.astype(COLUMN_TYPES['known_splice_junction'])
        all_combined_SJ_out_tab_df.loc[:, 'known_splice_junction'] = known_splice_junction_df

        # merged other columns across samples in the current batch
        all_combined_SJ_out_tab_df.loc[:, 'unique_reads'] = all_combined_SJ_out_tab_df[current_batch_columns['unique_reads']].sum(axis=1).astype(COLUMN_TYPES['unique_reads'])
        all_combined_SJ_out_tab_df.loc[:, 'multi_mapped_reads'] = all_combined_SJ_out_tab_df[current_batch_columns['multi_mapped_reads']].sum(axis=1).astype(COLUMN_TYPES['multi_mapped_reads'])
        all_combined_SJ_out_tab_df.loc[:, 'total_reads'] = all_combined_SJ_out_tab_df[current_batch_columns['total_reads']].sum(axis=1).astype(COLUMN_TYPES['total_reads'])
        all_combined_SJ_out_tab_df.loc[:, 'maximum_overhang'] = all_combined_SJ_out_tab_df[current_batch_columns['maximum_overhang']].max(axis=1).astype(COLUMN_TYPES['maximum_overhang'])
        all_combined_SJ_out_tab_df.loc[:, 'strand_counter'] = all_combined_SJ_out_tab_df[current_batch_columns['strand_counter']].sum(axis=1).astype(COLUMN_TYPES['strand_counter'])

        # compute derived columns
        all_combined_SJ_out_tab_df.loc[:, 'num_samples_with_this_junction'] = all_combined_SJ_out_tab_df[current_batch_columns['num_samples_with_this_junction']].sum(axis=1).astype(COLUMN_TYPES['num_samples_with_this_junction'])
        all_combined_SJ_out_tab_df.loc[:, 'num_samples_total'] = total_samples
        all_combined_SJ_out_tab_df.loc[:, 'max_per_sample_unique_reads'] = all_combined_SJ_out_tab_df[current_batch_columns['unique_reads']].max(axis=1).astype(COLUMN_TYPES['unique_reads'])
        all_combined_SJ_out_tab_df.loc[:, 'max_per_sample_total_reads'] = all_combined_SJ_out_tab_df[current_batch_columns['total_reads']].max(axis=1).astype(COLUMN_TYPES['total_reads'])

        if not args.discard_sample_id_column:
            all_combined_SJ_out_tab_df.loc[:, 'sample_id'] = all_combined_SJ_out_tab_df[current_batch_columns['sample_id']].apply(lambda row: ",".join(row[~row.isna()]), axis=1).astype('str')

        for column in COLUMN_NAMES:
            if batch_number > 0: current_batch_columns[column].remove(column)

            #logging.info(f"Dropping columns: {current_batch_columns[column]} from {all_combined_SJ_out_tab_df.columns}")
            all_combined_SJ_out_tab_df.drop(columns=current_batch_columns[column], inplace=True)

        print_memory_stats(f'after table {sample_i}')
        logging.info(f"Done processing batch {batch_number}")

    if args.normalize_read_counts:
        # undo the scaling that was applied to unique read counts with the expectation that they'd be summed (instead of taking the max value as is done for this column)
        all_combined_SJ_out_tab_df.loc[:, 'max_per_sample_unique_reads'] *= total_samples
        all_combined_SJ_out_tab_df.loc[:, 'max_per_sample_total_reads'] *= total_samples

    return all_combined_SJ_out_tab_df


def main():
    args = parse_args()

    all_combined_SJ_out_tab_df = process_all_batches(args)

    # set final strand value to 1 (eg. '+') or 2 (eg. '-') or 0 (eg. uknown) based on the setting in the majority of samples
    all_combined_SJ_out_tab_df.loc[:, 'strand'] = all_combined_SJ_out_tab_df['strand_counter'].apply(
        lambda s: 1 if s > 0 else (2 if s < 0 else 0)).astype('int8')

    logging.info(all_combined_SJ_out_tab_df.dtypes)
    logging.info("-----")
    #pd.set_option('display.max_columns', 10000)
    #logging.info(all_combined_SJ_out_tab_df.describe())

    # save the combined table
    if args.output_path:
        output_path = args.output_path
    else:
        output_path = f"combined.{len(args.paths)}_samples"
        if args.normalize_read_counts:
            output_path += ".normalized"
        output_path += ".SJ.out.tsv.gz"

    out = all_combined_SJ_out_tab_df.reset_index().astype(COLUMN_TYPES)

    out[["chrom", "start_1based", "end_1based"] + COLUMN_NAMES].to_csv(output_path, index=False, sep="\t")
    logging.info(f"Wrote out {output_path}")

    # save per-sample tables
    if args.save_individual_tables:
        if args.discard_sample_id_column:
            AGGREGATED_COLUMNS_TO_ADD_TO_INDIVIDUAL_TABLES.remove('sample_id')

        for path in args.paths:
            df = read_SJ_out_tab(path)
            df.drop(columns=AGGREGATED_COLUMNS_TO_ADD_TO_INDIVIDUAL_TABLES, inplace=True)
            df = df.join(all_combined_SJ_out_tab_df[AGGREGATED_COLUMNS_TO_ADD_TO_INDIVIDUAL_TABLES], how="left")
            out = df.reset_index()
            output_path = re.sub("(.SJ.out)?(.tsv|.tab)(.gz)?$", "", os.path.basename(path))
            if args.normalize_read_counts:
                output_path += ".normalized"
            output_path += ".SJ.out.tsv.gz"

            logging.info(f"Wrote out {output_path}")
            out[["chrom", "start_1based", "end_1based"] + COLUMN_NAMES].to_csv(output_path, index=False, sep="\t")


if __name__ == "__main__":
    main()

"""
#strand_conflicts_count = combined_ht.filter(hl.abs(combined_ht.strand_counter)/hl.float(combined_ht.num_samples_with_this_junction) < 0.1, keep=True).count()
#joined_table.to_csv(f"combined_using_pandas.{total_samples}_samples.SJ.out.tab", sep="\t", header=True, index=False)
#joined_table = pd.read_parquet(f"combined_using_pandas.{total_samples}_samples.SJ.out.parquet")

#total_junctions_count = combined_ht.count()
#strand_conflicts_count = combined_ht.filter(hl.abs(combined_ht.strand_counter)/hl.float(combined_ht.num_samples_with_this_junction) < 0.1, keep=True).count()

#combined_ht = combined_ht.annotate_globals(
#    combined_tables=args.paths,
#    n_combined_tables=total_samples)

#if strand_conflicts_count:
#    print(f"WARNING: Found {strand_conflicts_count} strand_conflicts out of {total_junctions_count} total_junctions")
"""
