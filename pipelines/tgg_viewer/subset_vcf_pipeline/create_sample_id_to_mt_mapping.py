"""
23439131  2017-12-07T17:56:41Z  gs://seqr-datasets/GRCh38/1kg/1kg.liftover.vcf.gz
239969348564  2017-07-04T22:33:18Z  gs://seqr-datasets/GRCh38/20170629_900Genomes_full_239969348564/900Genomes_full.vcf.gz
 453659893  2017-10-21T03:47:26Z  gs://seqr-datasets/GRCh38/ATGU_WGS-Jueppner/ATGU_WGS-Jueppner.vcf.gz
2203880710  2017-11-15T14:50:03Z  gs://seqr-datasets/GRCh38/CMG_Gazda/MacArthur_Muscle_PCRfree_WGS_v10_rehead.subset.vcf.gz
2203880710  2017-11-15T12:12:02Z  gs://seqr-datasets/GRCh38/CMG_Gazda_WGS/MacArthur_Muscle_PCRfree_WGS_v10_rehead.subset.vcf.gz
1266738360  2017-11-05T18:11:05Z  gs://seqr-datasets/GRCh38/CMG_Manton_Genomes/MacArthur_Muscle_PCRfree_WGS_v10.subset.vcf.gz
61900720496  2017-09-26T00:49:10Z  gs://seqr-datasets/GRCh38/GMKF_DSD/GMKF_DSD_final.vcf.gz
7142664891  2017-09-18T21:33:18Z  gs://seqr-datasets/GRCh38/GMKF_sex_disorders/sra-1.vcf.gz
1591485285  2017-11-05T18:14:05Z  gs://seqr-datasets/GRCh38/MYOSEQ_PCRFree_WGS_v5/MacArthur_Muscle_PCRfree_WGS_v10.subset.vcf.gz
 655626800  2018-02-26T03:28:31Z  gs://seqr-datasets/GRCh38/INMR_v9/INMR_v9.liftover.b38.vcf.gz

289429266127  2018-06-05T13:21:57Z  gs://seqr-datasets/GRCh38/CMG_RGP_Broad_MacArthur_RareDisease_WGS_Callsets/v2/CMG_RGP_Broad_MacArthur_RareDisease_WGS_Callsets.vcf.gz
312664082637  2018-05-24T14:06:23Z  gs://seqr-datasets/GRCh38/CMG_RGP_Broad_MacArthur_RareDisease_WGS_Callsets/v3/CMG_RGP_Broad_MacArthur_RareDisease_WGS_Callsets.vcf.gz

 163837989  2019-02-21T20:21:32Z  gs://seqr-datasets/GRCh38/RDG_WGS_Broad_Internal/v5/RDG_WGS_Broad_Internal/sharded-vcfs/RDG_WGS_Broad_Internal.filtered.0.vcf.gz
 170061229  2019-04-02T14:51:19Z  gs://seqr-datasets/GRCh38/RDG_WGS_Broad_Internal/v6/sharded_vcf/RDG_WGS_Broad_Internal_v5.filtered.0.vcf.gz
 210692572  2019-07-09T14:03:17Z  gs://seqr-datasets/GRCh38/RDG_WGS_Broad_Internal/v7/sharded_vcf/RDG_WGS_Broad_Internal.filtered.0.vcf.gz
 184414308  2019-09-16T14:08:35Z  gs://seqr-datasets/v02/GRCh38/RDG_WGS_Broad_Internal/v8/sharded_vcf/RDG_WGS_Broad_Internal.filtered.0.vcf.gz

1921549555  2019-06-11T15:38:14Z  gs://seqr-datasets/GRCh38/RDG_WES_Broad_Internal/v11/sharded_vcf/RDG_WES_Broad_Internal.filtered.0.vcf.gz
 143735863  2019-08-05T18:33:58Z  gs://seqr-datasets/GRCh38/RDG_WES_Broad_Internal/v12/sharded_vcf/RDG_Broad_WES_Internal.filtered.0.vcf.gz
226397595519  2019-08-27T14:27:14Z  gs://seqr-datasets/v02/GRCh38/RDG_WES_Broad_Internal/v13/RDG_Broad_WES_Internal.vcf.bgz
"""

#%%

import collections
import hail as hl
import os


grch38_vcfs = [
    ("gs://seqr-datasets/v02/GRCh38/RDG_WES_Broad_Internal/v16/RDG_WES_Broad_Internal.mt", False),

    ("gs://seqr-datasets/GRCh38/INMR_v9/INMR_v9.liftover.b38.vcf.gz", False),
    ("gs://seqr-datasets/GRCh38/MYOSEQ_PCRFree_WGS_v5/MacArthur_Muscle_PCRfree_WGS_v10.subset.vcf.gz", True),
    ("gs://seqr-datasets/GRCh38/CMG_Manton_Genomes/MacArthur_Muscle_PCRfree_WGS_v10.subset.vcf.gz", True),
    ("gs://seqr-datasets/GRCh38/CMG_Gazda/MacArthur_Muscle_PCRfree_WGS_v10_rehead.subset.vcf.gz", True),
    ("gs://seqr-datasets/GRCh38/20170629_900Genomes_full_239969348564/900Genomes_full.vcf.gz", True),
    ("gs://seqr-datasets/GRCh38/GMKF_DSD/GMKF_DSD_final.vcf.gz", False),

    ("gs://seqr-datasets/v02/GRCh38/RDG_WGS_Broad_Internal/v10/RDG_WGS_Broad_Internal.mt", False),
]

OUTPUT_DIR = "gs://seqr-bw/project__rnaseq/"

#%%

CONTIG_RECODING = {f"{c}": f"chr{c}" for c in list(range(1,23)) + ["X", "Y", "M"]}

sample_id_to_path = collections.OrderedDict()
for path, needs_vcf_header in grch38_vcfs:
    print("-----------")
    print(f"Importing {path}")
    dirname = os.path.dirname(path)
    if needs_vcf_header:
        mt = hl.import_vcf(path, reference_genome="GRCh38", force_bgz=True,
                           header_file=os.path.join(dirname, "header_for_hail_v0.2.vcf"),
                           contig_recoding=CONTIG_RECODING,
                           skip_invalid_loci=True)
    else:
        if path.endswith(".mt"):
            mt = hl.read_matrix_table(path)
        else:
            mt = hl.import_vcf(path, reference_genome="GRCh38", force_bgz=True, contig_recoding=CONTIG_RECODING, skip_invalid_loci=True)

    sample_ids = mt.cols().s.collect()
    print(f"Processing {len(sample_ids)} samples from {path}")
    counter = 0
    for sample_id in sample_ids:
        if sample_id not in sample_id_to_path:
            counter += 1
        sample_id_to_path[sample_id] = path
    print(f"Found {counter} new sample ids")


#%%

output_filename = "sample_id_to_vcf_path_mapping.tsv"
with open(output_filename, "wt") as f:
    f.write("sample_id\tvcf_path\n")
    counter = 0
    for sample_id, vcf_path in sorted(sample_id_to_path.items(), key=lambda x: x[0]):
        f.write(f"{sample_id}\t{vcf_path}\n")
        counter += 1

print(f"Wrote {counter} rows {output_filename}")


#%%

def subset_vcf(mt, dna_sample_id, vcf_output_path):
    mt = mt.filter_cols(mt.s == dna_sample_id)
    sample_ids = mt.cols().s.collect()
    if len(sample_ids) != 1:
        raise ValueError(f"ERROR: Couldn't find {dna_sample_id} in mt")

    mt = mt.filter_rows(hl.agg.any(mt.GT.is_non_ref()))

    hl.export_vcf(mt, vcf_output_path)

    return mt


#%%

joint_called_mt_path = VCF_TO_MT[joint_called_vcf_path]
print(f"Reading in {joint_called_mt_path}")

mt = hl.read_matrix_table(joint_called_mt_path)

SAMPLE_ID_MAP2 = {}
for rna_sample_id, dna_sample_id, vcf_path in SAMPLE_ID_MAP.values():
    single_sample_vcf_path = subset_vcf(vcf_path, rna_sample_id, dna_sample_id)
    SAMPLE_ID_MAP2[rna_sample_id] = (rna_sample_id, dna_sample_id, single_sample_vcf_path)


snps_output_path = f"gs://macarthurlab-rnaseq-data/grch38_vcfs/{rna_sample_id.replace(' ', '_')}.SNPs.vcf.gz"

print(f"{mt.count()} {rna_sample_id}  {dna_sample_id}   {snps_output_path}")


