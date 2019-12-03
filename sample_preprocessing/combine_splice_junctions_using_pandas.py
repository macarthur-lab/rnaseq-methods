import argparse
import gc
import logging
import pandas as pd
import numpy as np
import psutil
import os

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s: %(message)s')


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("-N", "--normalize-read-counts", action="store_true", help="whether to normalize read counts")
    p.add_argument("paths", nargs="+", help="Paths of 1 or more SJ.out.tab tables")
    return p.parse_args()


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


def main():
    args = parse_args()

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

    for i, path in enumerate(args.paths):

        """
        column 1: chromosome
        column 2: first base of the intron (1-based)
        column 3: last base of the intron (1-based)
        column 4: strand (0: undefined, 1: +, 2: -)
        column 5: intron motif: 0: non-canonical; 1: GT/AG, 2: CT/AC, 3: GC/AG, 4: CT/GC, 5:AT/AC, 6: GT/AT
        column 6: 0: unannotated, 1: annotated (only if splice junctions database is used)
        column 7: number of uniquely mapping reads crossing the junction
        column 8: number of multi-mapping reads crossing the junction
        column 9: maximum spliced alignment overhang
        """
        print(f"Table {i}: {path}")
        df = pd.read_table(
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
            dtype=COLUMN_TYPES)

        df['num_samples_with_this_junction'] = np.int32(1)
        df['strand_counter'] = df['strand'].apply(lambda s: np.int8(1 if s == 1 else (-1 if s == 2 else 0)))

        if i == 0:
            df['unique_reads_0'] = df['unique_reads']
            df['multi_mapped_reads_0'] = df['multi_mapped_reads']
            joined_table = df

            print_memory_stats(f'after table {i}')
            continue

        joined_table = joined_table.join(df, how="outer", rsuffix=f"_{i}")

        ## in rare cases, the strand for the same junction may differ across samples, so use a 2-step process that assigns strand based on majority of samples
        #joined_table.loc[joined_table['strand'].isnull(), 'strand'] = joined_table[f'strand_{i}']
        #joined_table['strand'] = joined_table['strand'].astype(COLUMN_TYPES['strand'])

        joined_table.loc[joined_table['intron_motif'].isnull(), f'intron_motif'] = joined_table[f'intron_motif_{i}']
        joined_table['intron_motif'] = joined_table['intron_motif'].astype(COLUMN_TYPES['intron_motif'])

        joined_table['known_splice_junction'] = (joined_table['known_splice_junction'].fillna(False).astype('int8') | joined_table[f'known_splice_junction_{i}'].fillna(False).astype('int8'))
        joined_table['unique_reads'] = joined_table[['unique_reads', f'unique_reads_{i}']].sum(axis=1).fillna(0).astype(COLUMN_TYPES['unique_reads'])
        joined_table['multi_mapped_reads'] = joined_table[['multi_mapped_reads', f'multi_mapped_reads_{i}']].sum(axis=1).fillna(0).astype(COLUMN_TYPES['multi_mapped_reads'])
        joined_table['maximum_overhang'] = joined_table[['maximum_overhang', f'maximum_overhang_{i}']].max(axis=1).fillna(0).astype(COLUMN_TYPES['maximum_overhang'])
        joined_table['num_samples_with_this_junction'] = joined_table[['num_samples_with_this_junction', f'num_samples_with_this_junction_{i}']].sum(axis=1).astype(COLUMN_TYPES['num_samples_with_this_junction'])
        joined_table['strand_counter'] = joined_table[['strand_counter', f'strand_counter_{i}']].sum(axis=1).fillna(0).astype(COLUMN_TYPES['strand_counter'])

        joined_table.drop(columns=[f'strand_{i}', f'strand_counter_{i}', f'intron_motif_{i}', f'known_splice_junction_{i}', f'maximum_overhang_{i}', f'num_samples_with_this_junction_{i}'], inplace=True)

        print_memory_stats(f'after table {i}')

    # set final strand value to 1 (eg. '+') or 2 (eg. '-') or 0 (eg. uknown) based on the setting in the majority of samples
    joined_table['strand'] = joined_table['strand_counter'].apply(lambda s: 1 if s > 0 else (2 if s < 0 else 0)).astype('int8')

    joined_table = joined_table.reset_index()

    read_count_columns = [f'unique_reads_{i}' for i in range(len(args.paths))] + [f'multi_mapped_reads_{i}' for i in range(len(args.paths))]
    joined_table[read_count_columns] = joined_table[read_count_columns].fillna(0).astype('int32')

    logging.info(joined_table.dtypes)
    logging.info("-----")
    pd.set_option('display.max_columns', 500)
    logging.info(joined_table.describe())

    joined_table.to_parquet(f"combined_using_pandas.{len(args.paths)}_samples.SJ.out.parquet")


if __name__ == "__main__":
    main()

"""
#strand_conflicts_count = combined_ht.filter(hl.abs(combined_ht.strand_counter)/hl.float(combined_ht.num_samples_with_this_junction) < 0.1, keep=True).count()
#joined_table.to_csv(f"combined_using_pandas.{len(args.paths)}_samples.SJ.out.tab", sep="\t", header=True, index=False)
#joined_table = pd.read_parquet(f"combined_using_pandas.{len(args.paths)}_samples.SJ.out.parquet")


ht = ht.annotate(
    strand_counter=hl.or_else(hl.switch(ht.strand).when(1, 1).when(2, -1).or_missing(), 0),
    num_samples_with_this_junction=1,
)

combined_ht = combined_ht.join(ht, how="outer")
combined_ht = combined_ht.transmute(
    strand=hl.or_else(combined_ht.strand, combined_ht.strand_1), ## in rare cases, the strand for the same junction may differ across samples, so use a 2-step process that assigns strand based on majority of samples
    strand_counter=hl.sum([combined_ht.strand_counter, combined_ht.strand_counter_1]),  # samples vote on whether strand = 1 (eg. '+') or 2 (eg. '-')
    intron_motif=hl.or_else(combined_ht.intron_motif, combined_ht.intron_motif_1),  ## double-check that left == right?

    known_splice_junction=hl.or_else(hl.cond((combined_ht.known_splice_junction == 1) | (combined_ht.known_splice_junction_1 == 1), 1, 0), 0), ## double-check that left == right?
    unique_reads=hl.sum([combined_ht.unique_reads, combined_ht.unique_reads_1]),
    multi_mapped_reads=hl.sum([combined_ht.multi_mapped_reads, combined_ht.multi_mapped_reads_1]),
    maximum_overhang=hl.max([combined_ht.maximum_overhang, combined_ht.maximum_overhang_1]),
    num_samples_with_this_junction=hl.sum([combined_ht.num_samples_with_this_junction, combined_ht.num_samples_with_this_junction_1]),
)

combined_ht = combined_ht.checkpoint(f"checkpoint{i % 2}.ht", overwrite=True) #, _read_if_exists=True)


#total_junctions_count = combined_ht.count()
#strand_conflicts_count = combined_ht.filter(hl.abs(combined_ht.strand_counter)/hl.float(combined_ht.num_samples_with_this_junction) < 0.1, keep=True).count()

# set final strand value to 1 (eg. '+') or 2 (eg. '-') or 0 (eg. uknown) based on the setting in the majority of samples
#combined_ht = combined_ht.annotate(
#    strand=hl.case()
#        .when(combined_ht.strand_counter > 0, 1)
#        .when(combined_ht.strand_counter < 0, 2)
#        .default(0))

#combined_ht = combined_ht.annotate_globals(
#    combined_tables=args.paths,
#    n_combined_tables=len(args.paths))

#if strand_conflicts_count:
#    print(f"WARNING: Found {strand_conflicts_count} strand_conflicts out of {total_junctions_count} total_junctions")

# write as HT
#combined_ht = combined_ht.checkpoint(f"combined.SJ.out.ht", overwrite=True) #, _read_if_exists=True)
"""
