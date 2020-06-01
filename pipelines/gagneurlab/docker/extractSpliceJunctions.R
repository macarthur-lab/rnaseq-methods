library(FRASER)
library(data.table)
library(stringr)
library(purrr)
library(argparse)

parser = ArgumentParser()
parser$add_argument("bam_header_path", help="Bam file for getting chromosome names")
args = parser$parse_args()

file_paths = list.files('.', pattern = "fraser_count_split_reads_.*.tar.gz$")
print(file_paths)
parse_sample_id = function(x) { return( str_replace(x[[1]], 'fraser_count_split_reads_', '')) }
sample_ids = unlist(map(strsplit(file_paths, '[.]'), parse_sample_id))

sampleTable = data.table(sampleID=sample_ids, bamFile=args$bam_header_path)
print(sampleTable)

fds = FraserDataSet(colData=sampleTable, workingDir=".", bamParam=ScanBamParam(mapqFilter=0), strandSpecific=0L)
splitCountsForAllSamples = getSplitReadCountsForAllSamples(fds)
splitCountRanges = rowRanges(splitCountsForAllSamples)
print(splitCountRanges)

saveRDS(splitCountRanges, "spliceJunctions.RDS")
