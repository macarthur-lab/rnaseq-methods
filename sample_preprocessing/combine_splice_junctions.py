import argparse
import hail as hl


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("-N", "--normalize-read-counts", action="store_true", help="whether to normalize read counts")
    p.add_argument("paths", nargs="+", help="Paths of 1 or more SJ.out.tab tables")
    return p.parse_args()


def import_SJ_out_tab(path):
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

    ht = hl.import_table(
        path,
        no_header=True,
        impute=True,
    ).rename({
        "f0": "chrom",
        "f1": "start_1based",
        "f2": "end_1based",
        "f3": "strand",
        "f4": "intron_motif",
        "f5": "known_splice_junction",
        "f6": "unique_reads",
        "f7": "multi_mapped_reads",
        "f8": "maximum_overhang",
    })
    
    return ht


def print_stats(path, ht):

    # print some stats
    total_splice_junctions = ht.count()
    known_splice_junctions = ht.filter(ht.known_splice_junction == 1).count()
    novel_splice_junctions = ht.filter(ht.known_splice_junction == 0).count()
    
    assert known_splice_junctions + novel_splice_junctions == total_splice_junctions

    print(f"Table: {path}")
    print(f"{total_splice_junctions} total splice junctions")
    print(f"{novel_splice_junctions} novel splice junctions ({100*novel_splice_junctions/total_splice_junctions:0.1f}%)")


def main():
    args = parse_args()

    tables = []
    for i, path in enumerate(args.paths):

        ht = import_SJ_out_tab(path)
        ht = ht.key_by("chrom", "start_1based", "end_1based")

        ht = ht.annotate_globals(
            path = path,
            unique_reads_in_sample = ht.aggregate(hl.agg.sum(ht.unique_reads)),
            multi_mapped_reads_in_sample = ht.aggregate(hl.agg.sum(ht.multi_mapped_reads)),
        )

        # add 'interval' column
        #ht = ht.annotate(interval=hl.interval(
        #    hl.locus(ht.chrom, ht.start_1based, reference_genome=reference_genome),
        #    hl.locus(ht.chrom, ht.end_1based, reference_genome=reference_genome),))

        tables.append(ht)

    # compute mean
    mean_unique_reads_in_sample = sum([hl.eval(ht.unique_reads_in_sample) for ht in tables])/float(len(tables))
    mean_multi_mapped_reads_in_sample = sum([hl.eval(ht.multi_mapped_reads_in_sample) for ht in tables])/float(len(tables))
    print(f"mean_unique_reads_in_sample: {mean_unique_reads_in_sample:01f}, mean_multi_mapped_reads_in_sample: {mean_multi_mapped_reads_in_sample:01f}")

    combined_ht = None
    for i, ht in enumerate(tables):
        unique_reads_multiplier = mean_unique_reads_in_sample / float(hl.eval(ht.unique_reads_in_sample))
        multi_mapped_reads_multiplier = mean_multi_mapped_reads_in_sample / float(hl.eval(ht.multi_mapped_reads_in_sample))
        print(f"unique_reads_multiplier: {unique_reads_multiplier:01f}, multi_mapped_reads_multiplier: {multi_mapped_reads_multiplier:01f}")
        ht = ht.annotate(
            num_samples_with_this_junction = 1,
        )

        if args.normalize_read_counts:
            ht = ht.annotate(
                unique_reads=ht.unique_reads * unique_reads_multiplier,
                multi_mapped_reads=ht.multi_mapped_reads * multi_mapped_reads_multiplier
            )

        if combined_ht is None:
            combined_ht = ht
            continue

        print("----")
        print_stats(path, ht)
        
        combined_ht = combined_ht.join(ht, how="outer")
        combined_ht = combined_ht.transmute(
            #strand=hl.or_else(combined_ht.strand, combined_ht.strand_1), ## in rare cases, the strand for the same junction may differ across samples, so use a 2-step process that assigns strand based on majority of samples
            plus_vs_minus_strand=hl.or_else(combined_ht.plus_vs_minus_strand, 0) + hl.or_else(hl.switch(combined_ht.strand).when(1, 1).when(2, -1).or_missing(), 0),  # samples vote on whether strand = 1 (eg. '+') or 2 (eg. '-')
            intron_motif=hl.or_else(combined_ht.intron_motif, combined_ht.intron_motif_1),  ## double-check that left == right?
            known_splice_junction=(hl.cond((combined_ht.known_splice_junction == 1) | (combined_ht.known_splice_junction_1 == 1), 1, 0)), ## double-check that left == right?
            unique_reads=hl.sum([combined_ht.unique_reads, combined_ht.unique_reads_1]),
            multi_mapped_reads=hl.sum([combined_ht.multi_mapped_reads, combined_ht.multi_mapped_reads_1]),
            maximum_overhang=hl.max([combined_ht.maximum_overhang, combined_ht.maximum_overhang_1]),
            num_samples_with_this_junction=hl.sum([combined_ht.num_samples_with_this_junction, combined_ht.num_samples_with_this_junction_1]),
        )
        #combined_ht = combined_ht.drop('unique_reads_in_sample_1', 'multi_mapped_reads_in_sample_1')

        combined_ht = combined_ht.checkpoint(f"checkpoint{i % 2}.ht", overwrite=True) #, _read_if_exists=True)
    
    total_junctions_count = combined_ht.count()
    strand_conflicts_count = combined_ht.filter(hl.abs(combined_ht.plus_vs_minus_strand)/hl.float(combined_ht.num_samples_with_this_junction) < 0.1, keep=True).count()

    # set final strand value to 1 (eg. '+') or 2 (eg. '-') or 0 (eg. uknown) based on the setting in the majority of samples
    combined_ht = combined_ht.annotate(
        strand=hl.case()
            .when(combined_ht.plus_vs_minus_strand > 0, 1)
            .when(combined_ht.plus_vs_minus_strand < 0, 2)
            .default(0))

    if strand_conflicts_count:
        print(f"WARNING: Found {strand_conflicts_count} strand_conflicts out of {total_junctions_count} total_junctions")

    # write as HT
    combined_ht = combined_ht.checkpoint(f"combined.SJ.out.ht", overwrite=True) #, _read_if_exists=True)

    ## write as tsv
    combined_ht = combined_ht.key_by()
    combined_ht.export("combined.SJ.out.with_header.tab", header=True)

    combined_ht = combined_ht.select(
        "chrom",
        "start_1based",
        "end_1based",
        "strand",
        "intron_motif",
        "known_splice_junction",
        "unique_reads",
        "multi_mapped_reads",
        "maximum_overhang",
    )
    combined_ht.export("combined.SJ.out.tab", header=False)

    print(f"unique_reads_in combined table: {combined_ht.aggregate(hl.agg.sum(combined_ht.unique_reads))}")

    ## use zip-join
    #combined_ht2 = hl.Table.multi_way_zip_join(tables, "data", "globals")
    #combined_ht2 = combined_ht2.annotate(
    #    strand=hl.agg.array_agg(lambda elem: hl.agg.max(elem.strand), combined_ht2.data))
    #combined_ht2.export("SJ.out.combined2.tab", header=False)
    #print(combined_ht2.describe())

if __name__ == "__main__":
    main()
