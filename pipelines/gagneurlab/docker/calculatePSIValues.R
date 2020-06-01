library(FRASER)
library(data.table)
library(stringr)
library(purrr)
library(argparse)

parser = ArgumentParser()
parser$add_argument("spliceJunctionsPath", help="RDS file containing splice junctions")
args = parser$parse_args()

splitCountRanges = readRDS(args$spliceJunctionsPath)
print(splitCountRanges)


file_paths = list.files('.', pattern = "fraser_count_split_reads_.*.tar.gz$")
print(file_paths)
parse_sample_id = function(x) { return( str_replace(x[[1]], 'fraser_count_split_reads_', '')) }
sample_ids = unlist(map(strsplit(file_paths, '[.]'), parse_sample_id))

sampleTable = data.table(sampleID=sample_ids, bamFile=args$bam_header_path)
print(sampleTable)

fds = FraserDataSet(colData=sampleTable, workingDir=".", bamParam=ScanBamParam(mapqFilter=0), strandSpecific=0L)

splitCountsForAllSamples = getSplitReadCountsForAllSamples(fds)
nonSplitCountsForAllSamples = getNonSplitReadCountsForAllSamples(fds, splitCountRanges)
fds = addCountsToFraserDataSet(fds, splitCountsForAllSamples, nonSplitCountsForAllSamples)
fds = calculatePSIValues(fds)

saveRDS(fds, "fdsWithPSIValues.RDS")
