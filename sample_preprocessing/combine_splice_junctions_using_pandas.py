import argparse
import gc
import itertools
import logging
import pandas as pd
import numpy as np
import psutil
import os

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
    'maximum_overhang': 'int16',
    'num_samples_with_this_junction': 'int32',
    'strand_counter': 'int32',
}


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("-b", "--batch-size", help="How many tables to merge at once. Bigger batch sizes use more memory.", type=int, default=100)
    p.add_argument("paths", nargs="+", help="Paths of 1 or more SJ.out.tab tables")
    return p.parse_args()


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


def read_table(path, i):
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

    df['num_samples_with_this_junction'] = np.int32(1)
    df['strand_counter'] = df['strand'].apply(lambda s: 1 if s == 1 else (-1 if s == 2 else 0)).astype('int32')

    name_mapping = dict(zip(df.columns, [f"{c}_{i}" for c in df.columns]))
    df.rename(columns=name_mapping, inplace=True)

    return df


def main():

    args = parse_args()

    result = pd.DataFrame()
    column_names = ['strand', 'intron_motif', 'known_splice_junction', 'unique_reads', 'multi_mapped_reads', 'maximum_overhang', 'num_samples_with_this_junction', 'strand_counter']
    i = 0
    for batch_number, batch in enumerate(batched_iter(args.paths, args.batch_size)):
        tables_in_batch = []
        batch_start_i = i
        for path in batch:
            logging.info(f"Batch {batch_number}, Table {i}: {path}")
            df = read_table(path, i)
            tables_in_batch.append(df)
            i += 1
        batch_end_i = i

        result = result.join(tables_in_batch, how="outer")

        batch_columns = {}
        for column in column_names:
            batch_columns[column] = [f'{column}_{k}' for k in range(batch_start_i, batch_end_i)]
            if batch_number > 0: batch_columns[column].append(column)
        result['intron_motif'] = result[batch_columns['intron_motif']].ffill(axis=1).iloc[:,-1].astype(COLUMN_TYPES['intron_motif'])
        result['known_splice_junction'] = result[batch_columns['known_splice_junction']].replace(0, np.nan).ffill(axis=1).iloc[:,-1].fillna(0).astype(COLUMN_TYPES['known_splice_junction'])

        result['unique_reads'] = result[batch_columns['unique_reads']].sum(axis=1).fillna(0).astype(COLUMN_TYPES['unique_reads'])
        result['multi_mapped_reads'] = result[batch_columns['multi_mapped_reads']].sum(axis=1).fillna(0).astype(COLUMN_TYPES['multi_mapped_reads'])
        result['maximum_overhang'] = result[batch_columns['maximum_overhang']].max(axis=1).fillna(0).astype(COLUMN_TYPES['maximum_overhang'])
        result['num_samples_with_this_junction'] = result[batch_columns['num_samples_with_this_junction']].sum(axis=1).astype(COLUMN_TYPES['num_samples_with_this_junction'])
        result['strand_counter'] = result[batch_columns['strand_counter']].sum(axis=1).fillna(0).astype(COLUMN_TYPES['strand_counter'])

        for column in column_names:
            if column in ['unique_reads', 'multi_mapped_reads']:
                continue  # keep the per-sample read count columns
            if batch_number > 0: batch_columns[column].remove(column)
            result.drop(columns=batch_columns[column], inplace=True)

        print_memory_stats(f'after table {i}')
        logging.info(result.dtypes)
        logging.info("-----")
        pd.set_option('display.max_columns', 30)
        logging.info(result.describe())

    # set final strand value to 1 (eg. '+') or 2 (eg. '-') or 0 (eg. uknown) based on the setting in the majority of samples
    result['strand'] = result['strand_counter'].apply(lambda s: 1 if s > 0 else (2 if s < 0 else 0)).astype('int8')

    result = result.reset_index()

    read_count_columns = [f'unique_reads_{i}' for i in range(len(args.paths))] + [f'multi_mapped_reads_{i}' for i in range(len(args.paths))]
    result[read_count_columns] = result[read_count_columns].fillna(0).astype('int32')

    logging.info(result.dtypes)
    logging.info("-----")
    pd.set_option('display.max_columns', 30)
    logging.info(result.describe())

    result.to_parquet(f"combined_using_pandas.{len(args.paths)}_samples.SJ.out.parquet")


if __name__ == "__main__":
    main()

"""
#strand_conflicts_count = combined_ht.filter(hl.abs(combined_ht.strand_counter)/hl.float(combined_ht.num_samples_with_this_junction) < 0.1, keep=True).count()
#joined_table.to_csv(f"combined_using_pandas.{len(args.paths)}_samples.SJ.out.tab", sep="\t", header=True, index=False)
#joined_table = pd.read_parquet(f"combined_using_pandas.{len(args.paths)}_samples.SJ.out.parquet")

#total_junctions_count = combined_ht.count()
#strand_conflicts_count = combined_ht.filter(hl.abs(combined_ht.strand_counter)/hl.float(combined_ht.num_samples_with_this_junction) < 0.1, keep=True).count()

#combined_ht = combined_ht.annotate_globals(
#    combined_tables=args.paths,
#    n_combined_tables=len(args.paths))

#if strand_conflicts_count:
#    print(f"WARNING: Found {strand_conflicts_count} strand_conflicts out of {total_junctions_count} total_junctions")
"""
