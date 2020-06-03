library(FRASER)
library(data.table)
library(stringr)
library(purrr)
library(argparse)
library(BiocParallel)

#x = sessionInfo()
#print(x)
#print(x$loadedOnly)
#exit()

parser = ArgumentParser()
#parser$add_argument("-n", "--num-threads", type="integer", default=1, help="Number of CPUs to use [default=%default]")
parser$add_argument("spliceJunctionsPath", help="RDS file containing splice junctions")
parser$add_argument("bamHeaderPath", help="Bam file for getting chromosome names")
args = parser$parse_args()

splitCountRanges = readRDS(args$spliceJunctionsPath)
print(splitCountRanges)


file_paths = list.files('.', pattern = "fraser_count_split_reads_.*.tar.gz$")
print(file_paths)
parse_sample_id = function(x) { return( str_replace(x[[1]], 'fraser_count_split_reads_', '')) }
sample_ids = unlist(map(strsplit(file_paths, '[.]'), parse_sample_id))

sampleTable = data.table(sampleID=sample_ids, bamFile=args$bamHeaderPath)
print(sampleTable)

fds = FraserDataSet(colData=sampleTable, workingDir=".", bamParam=ScanBamParam(mapqFilter=0), strandSpecific=0L)
bpparam = SerialParam(log=TRUE, progressbar=FALSE)
#bpparam = MulticoreParam(args$num_threads, log=TRUE, progressbar=FALSE)
splitCountsForAllSamples = getSplitReadCountsForAllSamples(fds, BPPARAM=bpparam)
nonSplitCountsForAllSamples = getNonSplitReadCountsForAllSamples(fds, splitCountRanges, BPPARAM=bpparam)
fds = addCountsToFraserDataSet(fds, splitCountsForAllSamples, nonSplitCountsForAllSamples)
fds = calculatePSIValues(fds)

saveRDS(fds, "fdsWithPSIValues.RDS")
