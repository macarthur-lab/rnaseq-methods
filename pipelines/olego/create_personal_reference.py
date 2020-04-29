import argparse
import hailtop.pipeline as hp
import os

BATCH_BILLING_PROJECT = "tgg-rare-disease"
GCLOUD_PROJECT = "seqr-project"
GCLOUD_USER_ACCOUNT = "weisburd@broadinstitute.org"
GCLOUD_CREDENTIALS_LOCATION = "gs://weisburd-misc/creds"


#HG38_FASTA = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"

def parse_args():
    p = argparse.ArgumentParser()
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument("--local", action="store_true", help="Batch: run locally")
    g.add_argument("--cluster", action="store_true", help="Batch: submit to cluster")
    p.add_argument("--name", help="Batch: (optional) job name")
    args = p.parse_args()

    return args


def main():
    args = parse_args()
    backend = hp.LocalBackend(gsa_key_file=os.path.abspath("misc-270914-cb9992ec9b25.json")) if args.local else hp.BatchBackend(args.project)
    p = hp.Pipeline(backend=backend, name=args.name)
    for i, vcf in enumerate([
        "gs://macarthurlab-rnaseq/grch38_vcfs/RGP_273_3_R1.SNPs.vcf.gz",
        "gs://macarthurlab-rnaseq/grch38_vcfs/RGP_248_3.SNPs.vcf.gz",
        "gs://macarthurlab-rnaseq/grch38_vcfs/RGP_54_3_2.SNPs.vcf.gz",

        #"gs://macarthurlab-rnaseq/grch38_vcfs/RGP_7_1_2.SNPs.vcf.gz",
        #"gs://macarthurlab-rnaseq/grch38_vcfs/RGP_7_2_2.SNPs.vcf.gz",
        #"gs://macarthurlab-rnaseq/grch38_vcfs/RGP_7_3_2.SNPs.vcf.gz",
        #"gs://macarthurlab-rnaseq/grch38_vcfs/RGP_7_4_2.SNPs.vcf.gz",
        #"gs://macarthurlab-rnaseq/grch38_vcfs/RGP_7_5_2.SNPs.vcf.gz",
    ]):

        t = p.new_task(name=f"olego{i}")

        hg38_fasta = p.read_input_group(fa=os.path.expanduser("~/p1/ref/GRCh38/hg38.fa"), fai=os.path.expanduser("~/p1/ref/GRCh38/hg38.fa.fai"))
        t.image("weisburd/olego:latest")

        # switch to user account
        t.command(f"gcloud auth activate-service-account --key-file /gsa-key/key.json")
        t.command(f"gsutil -m cp -r gs://weisburd-misc/creds/.config /tmp/")
        t.command(f"rm -rf ~/.config")
        t.command(f"mv /tmp/.config ~/")
        t.command(f"gcloud config set account weisburd@broadinstitute.org")
        t.command(f"gcloud config set project seqr-project")

        output_dir = os.path.basename(vcf).replace(".SNPs.vcf.gz", "")
        #t.command(f"gsutil -m cp {hg38_fasta} .")
        #t.command(f"gsutil -m cp {HG38_FASTA}.fai .")
        t.command(f"gsutil -m cp {vcf}* .")
        t.command(f"mkdir {output_dir}`")
        t.command(f"create_genomes --fasta {hg38_fasta.fa} --vcf {os.path.basename(vcf)} --outdir {output_dir}")
        t.command(f"gsutil -m cp -r {output_dir} gs://macarthurlab-rnaseq/grch38_personal_reference/")
        #if args.local:
        #    p.write_output(t.ofile, 'temp.txt')
        #else:
        #    p.write_output(t.ofile, 'gs://gnomad-bw2/temp.txt')

        p.run()
        break

    #t.command(f"gsutil ls > {t.ofile}")
    if isinstance(backend, hp.BatchBackend):
        backend.close()


if __name__ == "__main__":
    main()
