import hail as hl
import os

hl.init(log="/dev/null")

SAMPLES_TO_PROCESS = [
    ("./blood/combined.whole_blood.4_samples.SJ.out.tsv", "./blood/combined.gtex_blood.755_samples.SJ.out.tsv"),

    ("./lymph/combined.lymphocytes_M.2_samples.SJ.out.tsv", "./lymph/combined.gtex_lymphocytes.174_samples.SJ.out.tsv"),

    ("./fibs/combined.fibroblasts.84_samples.SJ.out.tsv", "./fibs/combined.gtex_fibs.504_samples.SJ.out.tsv"),
    ("./fibs/combined.fibroblasts_M.54_samples.SJ.out.tsv", "./fibs/combined.gtex_fibs.504_samples.SJ.out.tsv"),
    ("./fibs/combined.fibroblasts_F.30_samples.SJ.out.tsv", "./fibs/combined.gtex_fibs.504_samples.SJ.out.tsv"),

    ("./muscle/combined.muscle.154_samples.SJ.out.tsv", "./muscle/combined.gtex_muscle.803_samples.SJ.out.tsv"),
    ("./muscle/combined.muscle_F.51_samples.SJ.out.tsv", "./muscle/combined.gtex_muscle.803_samples.SJ.out.tsv"),
    ("./muscle/combined.muscle_M.103_samples.SJ.out.tsv", "./muscle/combined.gtex_muscle.803_samples.SJ.out.tsv"),

    ("./muscle/combined.muscle_M_101bp.58_samples.SJ.out.tsv", "./muscle/combined.gtex_muscle.803_samples.SJ.out.tsv"),
    ("./muscle/combined.muscle_F_101bp.28_samples.SJ.out.tsv", "./muscle/combined.gtex_muscle.803_samples.SJ.out.tsv"),
    ("./muscle/combined.muscle_M_76bp.45_samples.SJ.out.tsv", "./muscle/combined.gtex_muscle.803_samples.SJ.out.tsv"),
    ("./muscle/combined.muscle_F_76bp.23_samples.SJ.out.tsv", "./muscle/combined.gtex_muscle.803_samples.SJ.out.tsv"),



]

UPDATE_CHECKPOINTS = True


print(f"\n--------------------------------")
omim_path = "./OMIM_20201110.tsv.gz"
print(f"Reading {omim_path}")
omim_ht = hl.import_table(omim_path, force=True).repartition(200)
omim_ht = omim_ht.annotate(chrom="chr"+omim_ht.chrom)
print(f"Loaded {omim_ht.count()} OMIM records")

omim_ht = omim_ht.annotate(
    interval=hl.interval(
        start=hl.locus(omim_ht.chrom, hl.int(omim_ht.start), reference_genome="GRCh38"),
        end  =hl.locus(omim_ht.chrom, hl.int(omim_ht.end),   reference_genome="GRCh38"),
    )
).key_by("interval")

omim_ht = omim_ht.checkpoint("./checkpoints/omim_ht_checkpoint", _read_if_exists=not UPDATE_CHECKPOINTS, overwrite=True)

for input_samples_tsv, input_gtex_tsv in SAMPLES_TO_PROCESS:
    output_prefix =  os.path.join(os.path.dirname(input_samples_tsv), "_".join(os.path.basename(input_samples_tsv).split(".")[1:3]))
    print(f"\n--------------------------------")
    print(f"Processing {input_samples_tsv}")
    print(f"Output prefix: {output_prefix}")
    assert os.path.isfile(input_samples_tsv), input_samples_tsv
    assert os.path.isfile(input_gtex_tsv), input_gtex_tsv

    gtex_junctions_ht = hl.import_table(input_gtex_tsv).key_by("chrom", "start_1based", "end_1based").repartition(200)
    gtex_junctions_ht = gtex_junctions_ht.drop("strand", "known_splice_junction", "intron_motif", "strand_counter")
    columns = list(gtex_junctions_ht.row)
    gtex_junctions_ht = gtex_junctions_ht.rename({c: f"gtex_{c}" for c in columns})
    print(f"\nGot {gtex_junctions_ht.count()} rows from {input_gtex_tsv}")

    sample_junctions_ht = hl.import_table(input_samples_tsv).key_by("chrom", "start_1based", "end_1based").repartition(200)
    sample_junctions_ht.count()
    print(f"\nGot {sample_junctions_ht.count()} rows from {input_samples_tsv}")

    junctions_ht = sample_junctions_ht.join(gtex_junctions_ht, how="left").repartition(200)
    junctions_ht = junctions_ht.annotate(
        sample_id=junctions_ht.sample_id.split(","),
        junction_start_1based_locus=hl.locus(
            junctions_ht.chrom,
            hl.int(junctions_ht.start_1based),
            reference_genome="GRCh38"
        ),
        junction_end_1based_locus=hl.locus(
            junctions_ht.chrom,
            hl.int(junctions_ht.end_1based),
            reference_genome="GRCh38"
        ),
    ).key_by("junction_start_1based_locus")


    junctions_ht_not_in_omim = junctions_ht.filter(
        hl.is_defined(omim_ht[junctions_ht.junction_start_1based_locus]),
        keep=False,
    )

    junctions_ht_not_in_omim = junctions_ht_not_in_omim.annotate(
        omim=omim_ht[junctions_ht_not_in_omim.junction_start_1based_locus]
    ) # this is just so the schema matches junction_ht

    junctions_ht_not_in_omim = junctions_ht_not_in_omim.checkpoint(
        f"./checkpoints/{output_prefix}/junctions_not_in_omim_checkpoint.ht", _read_if_exists=not UPDATE_CHECKPOINTS, overwrite=True)
    junctions_ht_not_in_omim.count()


    junctions_ht_in_omim = junctions_ht.annotate(
        omim=omim_ht.index(junctions_ht.junction_start_1based_locus, all_matches=True)
    ).explode("omim")

    junctions_ht_in_omim = junctions_ht_in_omim.checkpoint(
        f"./checkpoints/{output_prefix}/junctions_in_omim_checkpoint.ht", _read_if_exists=not UPDATE_CHECKPOINTS, overwrite=True)
    junctions_ht.count()

    junctions_union_ht = junctions_ht_in_omim.union(junctions_ht_not_in_omim)
    junctions_union_ht = junctions_union_ht.checkpoint(
        "./checkpoints/{output_prefix}/junctions_union_checkpoint.ht", _read_if_exists=not UPDATE_CHECKPOINTS, overwrite=True)

    junctions_union_ht = junctions_union_ht.checkpoint(f"{output_prefix}.ht", _read_if_exists=not UPDATE_CHECKPOINTS, overwrite=True)

    # also save to .tsv
    junctions_union_ht = junctions_union_ht.annotate(
        gtex_unique_reads=hl.or_else(junctions_union_ht.gtex_unique_reads, "0"),
        gtex_multi_mapped_reads=hl.or_else(junctions_union_ht.gtex_multi_mapped_reads, "0"),
        gtex_maximum_overhang=hl.or_else(junctions_union_ht.gtex_maximum_overhang, "0"),
        gtex_num_samples_with_this_junction=hl.or_else(junctions_union_ht.gtex_num_samples_with_this_junction, "0"),
        sample_id=hl.delimit(junctions_union_ht.sample_id),
        junction_1based=junctions_union_ht.chrom + ":" + junctions_union_ht.start_1based + "-" + junctions_union_ht.end_1based,
    ).key_by().drop("junction_start_1based_locus", "junction_end_1based_locus")

    junctions_union_ht.describe()
    print(f"\n{output_prefix}  {junctions_union_ht.count()} rows")

    df = junctions_union_ht.to_pandas()
    df.to_csv(f"{output_prefix}.tsv", index=False, header=True, sep="\t")


