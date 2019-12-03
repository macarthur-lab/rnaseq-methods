import argparse
import gc
import pandas as pd
import psutil
import os


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("-N", "--normalize-read-counts", action="store_true", help="whether to normalize read counts")
    p.add_argument("paths", nargs="+", help="Paths of 1 or more SJ.out.tab tables")
    return p.parse_args()


def print_stats(path, ht):

    # print some stats
    total_splice_junctions = ht.count()
    known_splice_junctions = ht.filter(ht.known_splice_junction == 1).count()
    novel_splice_junctions = ht.filter(ht.known_splice_junction == 0).count()

    assert known_splice_junctions + novel_splice_junctions == total_splice_junctions

    print(f"Table: {path}")
    print(f"{total_splice_junctions} total splice junctions")
    print(f"{novel_splice_junctions} novel splice junctions ({100*novel_splice_junctions/total_splice_junctions:0.1f}%)")


def print_memory_stats(message="", run_gc=False):
    if run_gc:
        gc.collect()
    print(f'memory used - {message}: {psutil.Process(os.getpid()).memory_info().rss//10**6} Mb')


def main():
    args = parse_args()

    joined_table = None
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
            names=['chrom', 'start_1based', 'end_1based',
                    f'strand_{i}', f'intron_motif_{i}', f'known_splice_junction_{i}',
                    f'unique_reads_{i}', f'multi_mapped_reads_{i}', f'maximum_overhang_{i}'],
            index_col=['chrom', 'start_1based', 'end_1based'])

        if joined_table is None:
            joined_table = df
        else:
            joined_table = joined_table.join(df, how="outer")

        print_memory_stats(run_gc=True)

    joined_table = joined_table.reset_index()

    print_memory_stats('after reset index',  run_gc=True)

    #joined_table.to_csv(f"combined_using_pandas.{len(args.paths)}_samples.SJ.out.tab", sep="\t", header=True, index=False)
    joined_table.to_feather(f"combined_using_pandas.{len(args.paths)}_samples.SJ.out.feather")

    print_memory_stats('after exporting to feather', run_gc=True)

    joined_table.to_parquet(f"combined_using_pandas.{len(args.paths)}_samples.SJ.out.pqt")

    print_memory_stats('after exporting to parquet', run_gc=True)

    """
        ht = ht.annotate(
            strand_counter=hl.or_else(hl.switch(ht.strand).when(1, 1).when(2, -1).or_missing(), 0),
            num_samples_with_this_junction=1,
        )

        print("----")
        print_stats(path, ht)

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
    """

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


if __name__ == "__main__":
    main()
