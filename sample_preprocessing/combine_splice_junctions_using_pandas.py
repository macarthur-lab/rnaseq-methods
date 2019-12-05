import argparse
from collections import defaultdict
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
    'unique_reads': 'float32',
    'multi_mapped_reads': 'float32',
    'maximum_overhang': 'int16',
    'num_samples_with_this_junction': 'int32',
    'strand_counter': 'int32',
}


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("-b", "--batch-size", help="How many tables to merge at once. Bigger batch sizes use more memory.", type=int, default=100)
    p.add_argument("-N", "--normalize-read-counts", action="store_true", help="whether to normalize read counts")
    p.add_argument("-c", "--save-read-counts", action="store_true", help="Export separate parquet files with matrices of per-sample unique- and multi-mapped-read counts")
    p.add_argument("-o", "--output-path")
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


def read_table(path, i=None):
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

    # print some stats
    print(f"   {df.unique_reads.sum()/1_000_000:0.1f} million uniquely-mapped reads")
    print(f"   {100*df.multi_mapped_reads.sum()/float(df.unique_reads.sum() + df.multi_mapped_reads.sum()):0.0f}% multi-mapped")

    if i is not None:
        name_mapping = dict(zip(df.columns, [f"{c}_{i}" for c in df.columns]))
        df.rename(columns=name_mapping, inplace=True)

    return df


def write_to_parquet(df, path):
    df.to_parquet(path, index=True)
    logging.info(f"Wrote out {path}")


def main():

    args = parse_args()

    if args.normalize_read_counts:
        total_unique_reads_across_all_samples = 0
        for path in args.paths:
            df = read_table(path)
            total_unique_reads_across_all_samples += df.unique_reads.sum()

        average_unique_reads_per_sample = total_unique_reads_across_all_samples/float(len(args.paths))


    result = pd.DataFrame()
    column_names = ['strand', 'intron_motif', 'known_splice_junction', 'unique_reads', 'multi_mapped_reads', 'maximum_overhang', 'num_samples_with_this_junction', 'strand_counter']
    parquet_files = defaultdict(list)
    i = 0
    logging.info(f"Processing {len(args.paths)} tables")
    for batch_number, batch in enumerate(batched_iter(args.paths, args.batch_size)):
        tables_in_batch = []
        batch_start_i = i
        for path in batch:
            logging.info(f"Batch {batch_number}, Table {i}: {path}")
            df = read_table(path, i)
            if args.normalize_read_counts:
                unique_reads_in_sample = df[f"unique_reads_{i}"].sum()
                scalar = average_unique_reads_per_sample / float(unique_reads_in_sample)
                df[f"unique_reads_{i}"] *= scalar / float(len(args.paths))
                df[f"multi_mapped_reads_{i}"] *= scalar / float(len(args.paths))
                logging.info(f"{path} has {int(unique_reads_in_sample)} total unique reads while the per-sample average is {average_unique_reads_per_sample}. Scaling read counts by {scalar} and dividing by {len(args.path)}")
            tables_in_batch.append(df)
            i += 1
        batch_end_i = i

        result = result.join(tables_in_batch, how="outer")

        batch_columns = {}
        for column in column_names:
            batch_columns[column] = [f'{column}_{k}' for k in range(batch_start_i, batch_end_i)]
            if batch_number > 0: batch_columns[column].append(column)

        df = result[batch_columns['intron_motif']]
        df.ffill(axis=1, inplace=True)
        df = df.iloc[:,-1].astype(COLUMN_TYPES['intron_motif'])
        result['intron_motif'] = df

        df = result[batch_columns['known_splice_junction']]
        df.replace(0, np.nan, inplace=True)
        df.ffill(axis=1, inplace=True)
        df = df.iloc[:,-1]
        df.fillna(0, inplace=True)
        df = df.astype(COLUMN_TYPES['known_splice_junction'])
        result['known_splice_junction'] = df

        result['unique_reads'] = result[batch_columns['unique_reads']].sum(axis=1).astype(COLUMN_TYPES['unique_reads'])
        result['multi_mapped_reads'] = result[batch_columns['multi_mapped_reads']].sum(axis=1).astype(COLUMN_TYPES['multi_mapped_reads'])
        result['maximum_overhang'] = result[batch_columns['maximum_overhang']].max(axis=1).astype(COLUMN_TYPES['maximum_overhang'])
        result['num_samples_with_this_junction'] = result[batch_columns['num_samples_with_this_junction']].sum(axis=1).astype(COLUMN_TYPES['num_samples_with_this_junction'])
        result['strand_counter'] = result[batch_columns['strand_counter']].sum(axis=1).astype(COLUMN_TYPES['strand_counter'])

        for column in column_names:
            if batch_number > 0: batch_columns[column].remove(column)
            if args.save_read_counts and column in ['unique_reads', 'multi_mapped_reads']:
                output_file_name = f"{column}.batch_{batch_number}.{args.batch_size}_samples.SJ.out.parquet"
                parquet_files[column].append(output_file_name)

                read_count_df = result[batch_columns[column]].astype('float32')
                write_to_parquet(read_count_df, output_file_name)

            result.drop(columns=batch_columns[column], inplace=True)

        print_memory_stats(f'after table {i}')
        logging.info(f"Done processing batch {batch_number}")

    # set final strand value to 1 (eg. '+') or 2 (eg. '-') or 0 (eg. uknown) based on the setting in the majority of samples
    result['strand'] = result['strand_counter'].apply(lambda s: 1 if s > 0 else (2 if s < 0 else 0)).astype('int8')

    logging.info(result.dtypes)
    logging.info("-----")
    #pd.set_option('display.max_columns', 10000)
    #logging.info(result.describe())

    output_path = args.output_path or f"combined.{len(args.paths)}_samples"
    if args.normalize_read_counts:
        output_path += ".normalized"

    output_path += ".SJ.out.parquet"
    write_to_parquet(result, output_path)

    if args.save_read_counts:
        logging.info("Combine parquet files: ")
        for column in ['unique_reads', 'multi_mapped_reads']:
            tables = []
            for path in parquet_files[column]:
                tables.append(pd.read_parquet(path))
            result = pd.DataFrame().join(tables, how="outer")
            result.fillna(0, inplace=True)
            result = result.astype('float32')

            write_to_parquet(result, os.path.join(os.path.dirname(output_path), f"{column}.{os.path.basename(output_path)}"))


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
