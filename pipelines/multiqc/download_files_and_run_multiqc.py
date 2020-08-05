# Transfer files from a Terra bucket to the macarthurlab-rnaseq bucket

import argparse
import glob
import logging
import os
import pprint
import sys

logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s')
logger = logging.getLogger()
logger.setLevel(logging.INFO)


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--clean", action="store_true", help="Delete previous files")
    p.add_argument("batch_name", default="all", choices=["all", "batch_0", "batch_1_muntoni", "batch_2020_04", "batch_2020_08"])
    args = p.parse_args()

    return args


def run(cmd):
    logger.info(cmd)
    os.system(cmd)


def chdir(new_dir):
    logger.info("cd %s" % new_dir)
    os.chdir(new_dir)


def gsutil_cp(source, dest, mkdir=True):
    if mkdir:
        run("mkdir -p " + dest)
    run("gsutil -m cp -n %s  %s" % (source, dest))


def main():
    args = parse_args()
    logger.info("Args: " + pprint.pformat(args.__dict__))

    original_dir = os.path.abspath(os.getcwd())
    source_prefix = "gs://macarthurlab-rnaseq/" + (args.batch_name if args.batch_name != "all" else "*")
    dest_prefix = os.path.join(original_dir, "data/")
    if not os.path.isdir(dest_prefix):
        logger.error(dest_prefix + " dir not found.")
        sys.exit(1)

    dest_prefix = os.path.join(dest_prefix, args.batch_name)

    if args.clean:
        for subdirs in ["star", "rna_seqc", "fastqc"]:
            run("rm -r %s" % os.path.join(dest_prefix, subdirs))
        
    # star
    gsutil_cp("%s/star/*.Log.final.out" % source_prefix, dest_prefix+"/star/")
    gsutil_cp("%s/star/*.ReadsPerGene.out.tab.gz" % source_prefix,  dest_prefix+"/star/genecounts/")

    # rnaseqc
    gsutil_cp("%s/rnaseqc/*.metrics.tsv" % source_prefix,  dest_prefix+"/rna_seqc/metrics/")

    # fastqc
    gsutil_cp("%s/fastqc/zip/*_fastqc.zip" % source_prefix,  dest_prefix+"/fastqc/zip/")

    # unzip
    run("gunzip -f " + dest_prefix+"/star/genecounts/*.gz")
    chdir(dest_prefix+"/fastqc/zip/")
    for zip_path in glob.glob("*.zip"):
        run("unzip -o " + zip_path)

    # run multiqc
    chdir(dest_prefix)
    run("python3 -m multiqc -f -m star -m fastqc -m rna_seqc --filename %s.html ." % args.batch_name)

    run("rm -r %s" % os.path.join(original_dir, "%s_data" % args.batch_name))
    run("mv -f %s* %s" % (args.batch_name, original_dir))
    logger.info("Done")


if __name__ == "__main__":
    main()
