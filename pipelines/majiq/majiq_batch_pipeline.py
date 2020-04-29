import argparse
import hailtop.batch as hp
import os

from sample_metadata.utils import get_joined_metadata_df

GCLOUD_PROJECT = "seqr-project"
GCLOUD_USER_ACCOUNT = "weisburd@broadinstitute.org"
GCLOUD_CREDENTIALS_LOCATION = "gs://weisburd-misc/creds"


def parse_args():
    p = argparse.ArgumentParser()
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument("--local", action="store_true", help="Batch: run locally")
    g.add_argument("--cluster", action="store_true", help="Batch: submit to cluster")
    p.add_argument("--project", default="tgg-rare-disease", help="Batch: billing project. Required if submitting to cluster.")
    p.add_argument("--name", help="Batch: (optional) job name")
    args = p.parse_args()

    return args


def main():
    args = parse_args()

    df = get_joined_metadata_df()
    print("\n".join(df.columns))
    #print(", ".join(df.index))
    sample_ids = df[df['star_pipeline_batch'] == 'batch_1_muntoni'].sample_id

    print(f"Processing sample id sets: {', '.join(sample_ids)}")
    return

    # see https://hail.zulipchat.com/#narrow/stream/223457-Batch-support/topic/auth.20as.20user.20account for more details
    backend = hp.LocalBackend(gsa_key_file=os.path.expanduser("~/.config/gcloud/misc-270914-cb9992ec9b25.json")) if args.local else hp.ServiceBackend(args.project)
    p = hp.Batch(backend=backend, name=args.name)

    for sample_id in sample_ids:
        metadata_row = df.loc[sample_id]
        batch_name = metadata_row['star_pipeline_batch']

        t = p.new_job(name=args.name)
        t.image("weisburd/majiq:latest")

        if args.local:
            genes_gtf = p.read_input("gencode.v26.annotation.gff3", extension=".gff3")
            input_file = p.read_input_group(
                bam="250DV_LR_M1__100000_reads.Aligned.sortedByCoord.out.bam",  # metadata_row['star_bam'],
                bai="250DV_LR_M1__100000_reads.Aligned.sortedByCoord.out.bam.bai",  # metadata_row['star_bai'],
            )
        else:
            genes_gtf = p.read_input("gs://macarthurlab-rnaseq/ref/gencode.v26.annotation.GRCh38.gff3", extension=".gff3")
            input_file = p.read_input_group(
                bam=metadata_row['star_bam'],
                bai=metadata_row['star_bai'],
            )

        # switch to user account
        t.command(f"gcloud auth activate-service-account --key-file /gsa-key/key.json")
        t.command(f"gsutil -m cp -r {GCLOUD_CREDENTIALS_LOCATION}/.config /tmp/")
        t.command(f"rm -rf ~/.config")
        t.command(f"mv /tmp/.config ~/")
        t.command(f"gcloud config set account {GCLOUD_USER_ACCOUNT}")
        t.command(f"gcloud config set project {GCLOUD_PROJECT}")

        # run majiq build
        #t.command(f"gsutil -u {GCLOUD_PROJECT} -m cp gs://gtex-resources/GENCODE/gencode.v26.GRCh38.ERCC.genes.collapsed_only.gtf .")
        t.command(f"mv {genes_gtf} gencode.gff3")
        t.command(f"mv {input_file.bam} {sample_id}.bam")
        t.command(f"mv {input_file.bai} {sample_id}.bam.bai")

        t.command(f"echo '[info]' >> majiq_build.cfg")
        t.command(f"echo 'readlen={metadata_row['read length (rnaseqc)']}' >> majiq_build.cfg")
        t.command(f"echo 'bamdirs=.' >> majiq_build.cfg")
        t.command(f"echo 'genome=hg38' >> majiq_build.cfg")
        t.command(f"echo 'strandness={'None' if metadata_row['stranded? (rnaseqc)'] == 'no' else 'reverse'}' >> majiq_build.cfg")
        t.command(f"echo '[experiments]' >> majiq_build.cfg")
        t.command(f"echo '{sample_id}={sample_id}' >> majiq_build.cfg")

        t.command(f"cat majiq_build.cfg >> {t.ofile}")
        t.command(f"majiq build gencode.gff3 -c majiq_build.cfg -j 1 -o majiq_build_{sample_id} >> {t.ofile}")

        t.command(f"tar czf majiq_build_{sample_id}.tar.gz majiq_build_{sample_id}")
        t.command(f"cp majiq_build_{sample_id}.tar.gz {t.output_tar_gz}")

        t.command(f"ls -lh . >> {t.ofile}")
        t.command(f"echo ls majiq_build_{sample_id} >> {t.ofile}")
        t.command(f"ls -1 majiq_build_{sample_id} >> {t.ofile}")
        t.command(f"echo --- done gs://macarthurlab-rnaseq/{batch_name}/majiq_build/majiq_build_{sample_id}.tar.gz >> {t.ofile}")

        if args.local:
            p.write_output(t.output_tar_gz, f"majiq_build_{sample_id}.tar.gz")
            p.write_output(t.ofile, 'temp.txt')
        else:
            p.write_output(t.output_tar_gz, f"gs://macarthurlab-rnaseq/{batch_name}/majiq_build/majiq_build_{sample_id}.tar.gz")
            p.write_output(t.ofile, 'gs://gnomad-bw2/temp.txt')

    p.run()

    if isinstance(backend, hp.ServiceBackend):
        backend.close()


if __name__ == "__main__":
    main()


"""
Metadata columns:
star_pipeline_batch_x
star_bam
star_bai
star_SJ_out_tab
star_reads_per_gene_tab
grch38_vcf
rnaseqc_gene_reads
rnaseqc_exon_reads
rnaseqc_gene_tpm
rnaseqc_metrics
junctions_bed
coverage_bigwig
hg19_bam
hg19_bai
fastqc_zip
star_pipeline_batch_y
batch_date_from_hg19_bam_header
stranded? (rnaseqc)
read length (rnaseqc)
total reads x 10^6 (rnaseqc)
mapping rate (rnaseqc)
solved using RNA-seq?
solved not using RNA-seq?
proj (seqr)
fam (seqr)
proj2 (seqr)
fam2 (seqr)
sex
genome (seqr)
population (seqr)
sample type (seqr)
analysis status (seqr)
variant tags (seqr)
variant notes (seqr)
coded phenotype (seqr)
anlaysis summary + notes (seqr)
internal case review notes (seqr)
cram path (seqr)
Alias (Beryl)
Clinical Diagnosis (Beryl)
Sex (Beryl)
Age at muscle biopsy (Beryl)
Site of biopsy (Beryl)
Previous NGS testing (Beryl)
Genetic diagnosis before RNA-seq (Beryl)
Notes (Beryl)
"""
