import argparse
import os
import subprocess
import sys

if sys.platform == "darwin":
    platform = "macOSX"
elif sys.platform.startswith("linux"):
    platform = "linux"
else:
    platform = "[your OS]"

p = argparse.ArgumentParser(description="This script takes 2 or more bigWig files and merges them using bigWigMerge and bedGraphToBigWig. These tools can be downloaded from UCSC: " 
        f"http://hgdownload.soe.ucsc.edu/admin/exe/{platform}.x86_64/bigWigMerge     "
        f"http://hgdownload.soe.ucsc.edu/admin/exe/{platform}.x86_64/bigWigToBedGraph")

p.add_argument("-R", "--reference-fasta", default="hg38.fa.chromSizes", help="2-column text file created from the .fasta.fai by running 'cut -f 1,2 hg38.fa.fai > hg38.fa.chromSizes'")
p.add_argument('-o', '--output-path', help="Output .bigWig file path", default="merged.bigWig")
p.add_argument('input_path', help="Input .bigWig file path", nargs="+")
args = p.parse_args()

input_files = " ".join(args.input_path)

def run(c):
    print(c)
    subprocess.check_output(c, shell=True)


run(f"bigWigMerge {input_files} temp.bedGraph")
run("LC_COLLATE=C sort -k1,1 -k2,2n ./temp.bedGraph > temp2.bedGraph")
run("rm temp.bedGraph")
run(f"bedGraphToBigWig temp2.bedGraph {args.reference_fasta} {args.output_path}")
run("rm temp2.bedGraph")


