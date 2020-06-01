# based on docs @ https://bioconductor.org/packages/devel/bioc/vignettes/FRASER/inst/doc/FRASER.pdf

library(argparse)

parser = ArgumentParser()
parser$add_argument("sample_id", help="Sample ID")
parser$add_argument("bam_path", help="Bam file path")
parser$add_argument("splice_junctions_RDS_path", help="Path of RDS file containing rowRanges(splitCountsForAllSamples)")
args = parser$parse_args()

library(FRASER)
library(data.table)

spliceJunctions = readRDS(args$splice_junctions_RDS_path)

sampleTable = data.table(sampleID=c(args$sample_id), bamFile=c(args$bam_path))
fds = FraserDataSet(colData=sampleTable, workingDir=".", bamParam=ScanBamParam(mapqFilter=0), strandSpecific=0L)

getNonSplitReadCountsForAllSamples(fds, spliceJunctions)  # saves results to cache/



