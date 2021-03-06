# Transfer files from a Terra bucket to the tgg-rnaseq bucket

import argparse
import glob
import logging
import os
import pprint
import sys
from sample_metadata.rnaseq_metadata_utils import get_rnaseq_metadata_joined_with_paths_df

logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s')
logger = logging.getLogger()
logger.setLevel(logging.INFO)


def run(cmd):
    logger.info(cmd)
    os.system(cmd)


def chdir(new_dir):
    logger.info("cd %s" % new_dir)
    os.chdir(new_dir)


def gsutil_cp(source, dest, mkdir=True):
    if mkdir:
        run("mkdir -p " + dest)
    run("gsutil cp -n %s  %s" % (source, dest))


def main():
    rnaseq_sample_metadata_df = get_rnaseq_metadata_joined_with_paths_df()

    analysis_batches = set([b for b in rnaseq_sample_metadata_df["analysis batch"] if b.strip()])
    star_pipeline_batches = set([b for b in rnaseq_sample_metadata_df["star_pipeline_batch"] if b])

    p = argparse.ArgumentParser()
    p.add_argument("--clean", action="store_true", help="Delete previous files")
    p.add_argument("--dont-download", action="store_true", help="Skip the download files step")
    p.add_argument("batch_name", nargs="+", choices={"all",} | analysis_batches | star_pipeline_batches)
    args = p.parse_args()

    logger.info("Args: " + pprint.pformat(args.__dict__))

    for batch_name in args.batch_name:
        original_dir = os.path.abspath(os.getcwd())
        dest_prefix = os.path.join(original_dir, "data/")
        if not os.path.isdir(dest_prefix):
            logger.error(dest_prefix + " dir not found.")
            sys.exit(1)

        dest_prefix = os.path.join(dest_prefix, batch_name)
        run(f"mkdir -p {dest_prefix}")

        if batch_name == "all" or batch_name in star_pipeline_batches:
            if batch_name == "all":
                source_prefix = "gs://tgg-rnaseq/*"
            else:
                source_prefix = f"gs://tgg-rnaseq/{batch_name}"

            if not args.dont_download:
                # star
                gsutil_cp("%s/star/*.Log.final.out" % source_prefix, dest_prefix+"/star/")
                gsutil_cp("%s/star/*.ReadsPerGene.out.tab.gz" % source_prefix,  dest_prefix+"/star/genecounts/")

                # rnaseqc
                gsutil_cp("%s/rnaseqc/*.metrics.tsv" % source_prefix,  dest_prefix+"/rna_seqc/metrics/")

                # fastqc
                gsutil_cp("%s/fastqc/zip/*_fastqc.zip" % source_prefix,  dest_prefix+"/fastqc/zip/")

        elif batch_name in analysis_batches:
            df = rnaseq_sample_metadata_df[rnaseq_sample_metadata_df["analysis batch"] == batch_name]
            source_prefix = "gs://tgg-rnaseq/*"

            if not args.dont_download:
                star_logs = []
                star_gene_counts = []
                rnaseqc_metrics = []
                fastqc_zip = []
                for _, row in df.iterrows():
                    # star
                    star_logs.append(f"gs://tgg-rnaseq/{row.star_pipeline_batch}/star/{row.sample_id}*.Log.final.out")
                    star_gene_counts.append(f"gs://tgg-rnaseq/{row.star_pipeline_batch}/star/{row.sample_id}*.ReadsPerGene.out.tab.gz")

                    # rnaseqc
                    rnaseqc_metrics.append(f"gs://tgg-rnaseq/{row.star_pipeline_batch}/rnaseqc/{row.sample_id}*.metrics.tsv")

                    # fastqc
                    fastqc_zip.append(f"gs://tgg-rnaseq/{row.star_pipeline_batch}/fastqc/zip/{row.sample_id}*_fastqc.zip")

                gsutil_cp(" ".join(star_logs), dest_prefix+"/star/")
                gsutil_cp(" ".join(star_gene_counts), dest_prefix+"/star/genecounts/")
                gsutil_cp(" ".join(rnaseqc_metrics), dest_prefix+"/rna_seqc/metrics/")
                gsutil_cp(" ".join(fastqc_zip), dest_prefix+"/fastqc/zip/")
        else:
            raise ValueError(batch_name)

        #if args.clean:
        #    for subdirs in ["star", "rna_seqc", "fastqc"]:
        #        run("rm -r %s" % os.path.join(dest_prefix, subdirs))

        if not args.dont_download:
            # unzip
            run("gunzip -f " + dest_prefix+"/star/genecounts/*.gz")
            chdir(dest_prefix+"/fastqc/zip/")
            for zip_path in glob.glob("*.zip"):
                run("unzip -o " + zip_path)

        # run multiqc
        chdir(dest_prefix)
        run("python3 -m multiqc -f -m star -m fastqc -m rna_seqc --filename %s.html ." % batch_name)

        run("rm -r %s" % os.path.join(original_dir, "%s_data" % batch_name))
        run("mv -f %s* %s" % (batch_name, original_dir))
        logger.info("Done")

        chdir(original_dir)


if __name__ == "__main__":
    main()
