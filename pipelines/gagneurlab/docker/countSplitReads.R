# based on docs @ https://bioconductor.org/packages/devel/bioc/vignettes/FRASER/inst/doc/FRASER.pdf

library(argparse)

parser = ArgumentParser()
parser$add_argument("sample_id", help="Sample ID")
parser$add_argument("bam_path", help="Bam file path")
args = parser$parse_args()

library(FRASER)
library(data.table)

sampleTable = data.table(sampleID=c(args$sample_id), bamFile=c(args$bam_path))
fds = FraserDataSet(colData=sampleTable, workingDir=".", bamParam=ScanBamParam(mapqFilter=0), strandSpecific=0L)

getSplitReadCountsForAllSamples(fds)  # saves results to cache/
