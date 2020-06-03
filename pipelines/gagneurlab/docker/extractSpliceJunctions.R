library(FRASER)
library(data.table)
library(stringr)
library(purrr)
library(argparse)
library(BiocParallel)

parser = ArgumentParser()
#parser$add_argument("-n", "--num-threads", type="integer", default=1, help="Number of CPUs to use [default=%default]")
parser$add_argument("bamHeaderPath", help="Bam file for getting chromosome names")
args = parser$parse_args()

file_paths = list.files('.', pattern = "fraser_count_split_reads_.*.tar.gz$")
print(file_paths)
parse_sample_id = function(x) { return( str_replace(x[[1]], 'fraser_count_split_reads_', '')) }
sample_ids = unlist(map(strsplit(file_paths, '[.]'), parse_sample_id))

sampleTable = data.table(sampleID=sample_ids, bamFile=args$bamHeaderPath)
print(sampleTable)

bpparam = SerialParam(log=TRUE, progressbar=FALSE)
#bpparam = MulticoreParam(args$num_threads, log=TRUE, progressbar=FALSE)

fds = FraserDataSet(colData=sampleTable, workingDir=".", bamParam=ScanBamParam(mapqFilter=0), strandSpecific=0L)
splitCountsForAllSamples = getSplitReadCountsForAllSamples(fds, BPPARAM=bpparam)
splitCountRanges = rowRanges(splitCountsForAllSamples)
print(splitCountRanges)

saveRDS(splitCountRanges, "spliceJunctions.RDS")
