# Transfer files from a Terra bucket to the macarthurlab-rnaseq bucket

import argparse
import logging
import os
import pprint

from firecloud import api

logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s')
logger = logging.getLogger()
logger.setLevel(logging.INFO)

FILE_TYPES = ["hg19_bams", "star", "rnaseqc", "fastqc", "coverage"]

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--source-namespace", help="Terra namespace", default="macarthurlab-rnaseq-terra")
    p.add_argument("--source-workspace", help="Terra workspace name", default="macarthurlab-rnaseq-terra")
    p.add_argument("--dest-bucket", default="macarthurlab-rnaseq")
    p.add_argument("-w", "--workflow-id", help="(optional) workflow id. Can specify more than one", action="append")
    p.add_argument("-t", "--file-type", choices=FILE_TYPES, help="(optional) what types of files to transfer", action="append")
    p.add_argument("-f", "--force", action="store_true", help="Force copy even if destination files already exist.")
    p.add_argument("-d", "--delete-after-copy", action="store_true", help="Delete source files")

    p.add_argument("batch_name", choices=["batch_0", "batch_1_muntoni", "batch_2020_04", "batch_2020_08", "batch_2020_08__walsh", "batch_2021_01"])
    args = p.parse_args()

    if not args.file_type:
        args.file_type = FILE_TYPES

    return args


def run(cmd):
    logger.info(cmd)
    os.system(cmd)


def gsutil_cp(source, dest, force=False, delete_after_copy=False):
    n_arg = "" if force else "-n"
    delete_command = f" && gsutil -m rm -rf {source}" if delete_after_copy else ""

    run(f"gsutil cp {n_arg} {source} {dest} {delete_command}")


def copy_hg19_bams(args):
    # returns a list of dictionaries. Each dictionary is a row in one of the metadata tables.
    entities = api.get_entities_with_type(args.source_namespace, args.source_workspace).json()

    # if the workspace has both sample, sample_set, and participant tables, the list includes rows from all of them, so filter to just the sample records
    sample_entities = [e for e in entities if e['entityType'] == "sample"]
    # for each sample, get the relevant columns including cram path, and add them to a tsv file string for upload
    ws_total_bams = ws_total_bais = 0
    for e in sample_entities:
        sample_name = e['name']
        attr = e['attributes']
        # find which column in this workspace has the cram/bam path. Then copy it to the 'cram_path'/'crai_path' columns. I previously got this list of column names by retrieving all tables from all workspaces and looking for any column name with "cram" or "bam"
        found_bam_path = False
        for bam_key, bai_key in [
            ('bam_file', 'bai_file'),
            ('cram_or_bam_path', 'crai_or_bai_path'),
        ]:
            if attr.get(bam_key): #  and os.system("gsutil -u seqr-project ls " + attr[bam_key]) == 0:
                found_bam_path = True
                attr['hg19_bam_path'] = attr.get(bam_key)
                attr['hg19_bai_path'] = attr.get(bai_key)
                if attr.get('hg19_bam_path'):
                    break

        else:
            print("Missing bam path for: "  + sample_name)
            continue

        ws_total_bams += 1
        ws_total_bais += 1

        dest = "gs://%s/%s/hg19_bams/" % (args.dest_bucket, args.batch_name)
        gsutil_cp(attr['hg19_bam_path'], dest, delete_after_copy=args.delete_after_copy)
        gsutil_cp(attr['hg19_bam_path'].replace(".bam", "*bai"), dest, delete_after_copy=args.delete_after_copy)


def main():
    args = parse_args()
    logger.info("Args:\n" + pprint.pformat(args.__dict__))

    try:
        response = api.get_workspace(args.source_namespace, args.source_workspace)
    except Exception as e:
        logger.error(e)
        return

    try:
        ws = response.json()
    except Exception as e:
        logger.error(e)
        logger.error(response)
        return

    logger.info("Workspace:\n" + pprint.pformat(ws))

    source_bucket = ws['workspace']['bucketName']
    bucket_and_batch_name = (args.dest_bucket, args.batch_name)

    # hg19 bams
    #copy_hg19_bams(args)

    for workflow_id in (args.workflow_id or ['']):
        if workflow_id:
            logger.info("------------------------------------------------")
            logger.info("Processing workflow id: " + workflow_id)

        source_prefix = "gs://%s/%s" % (source_bucket, workflow_id)

        # star
        if "star" in args.file_type:
            dest = "gs://%s/%s/star/" % bucket_and_batch_name
            gsutil_cp("%s**star_out/*.Aligned.sortedByCoord.out.bam" % source_prefix, dest, force=args.force)
            gsutil_cp("%s**star_out/*.Aligned.sortedByCoord.out.bam.bai" % source_prefix,  dest, force=args.force)
            gsutil_cp("%s**star_out/*.Chimeric.out.junction.gz" % source_prefix, dest, force=args.force)
            gsutil_cp("%s**star_out/*.Log.final.out" % source_prefix,  dest, force=args.force)
            gsutil_cp("%s**star_out/*.Log.out" % source_prefix,  dest, force=args.force)
            gsutil_cp("%s**star_out/*.Log.progress.out" % source_prefix,  dest, force=args.force)
            gsutil_cp("%s**star_out/*.ReadsPerGene.out.tab.gz" % source_prefix,  dest, force=args.force)
            gsutil_cp("%s**star_out/*.SJ.out.tab.gz" % source_prefix,  dest, force=args.force)

        # rnaseqc
        if "rnaseqc" in args.file_type:
            dest = "gs://%s/%s/rnaseqc/" % bucket_and_batch_name
            gsutil_cp("%s**call-rnaseqc2/*.metrics.tsv" % source_prefix, dest, force=args.force)
            gsutil_cp("%s**call-rnaseqc2/*.exon_reads.gct.gz" % source_prefix, dest, force=args.force)
            gsutil_cp("%s**call-rnaseqc2/*.gene_reads.gct.gz" % source_prefix, dest, force=args.force)
            gsutil_cp("%s**call-rnaseqc2/*.gene_tpm.gct.gz" % source_prefix, dest, force=args.force)

        # fastqc
        if "fastqc" in args.file_type:
            dest = "gs://%s/%s/fastqc/zip/" % bucket_and_batch_name
            gsutil_cp("%s**_fastqc.zip" % source_prefix, dest, force=args.force)

        if "coverage" in args.file_type:
            dest = "gs://%s/%s/bigWig/" % bucket_and_batch_name
            gsutil_cp("%s**.bigWig" % source_prefix, dest, force=args.force)

    logger.info("Done")


if __name__ == "__main__":
    main()
