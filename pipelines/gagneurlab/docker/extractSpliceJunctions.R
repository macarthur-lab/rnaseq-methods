library(FRASER)
library(data.table)
library(stringr)
library(purrr)

file_paths = list.files('.', pattern = "fraser_count_rna.*.tar.gz$")
print(file_paths)
parse_sample_id = function(x) { return( str_replace(x[[1]], 'fraser_count_rna_', '')) }
sample_ids = unlist(map(strsplit(file_paths, '[.]'), parse_sample_id))

sampleTable = data.table(sampleID=sample_ids, bamFile=".")
print(sampleTable)

fds = FraserDataSet(colData=sampleTable, workingDir=".", bamParam=ScanBamParam(mapqFilter=0), strandSpecific=0L)
splitCountsForAllSamples = getSplitReadCountsForAllSamples(fds)
splitCountRanges = rowRanges(splitCountsForAllSamples)
saveRDS(splitCountRanges, "spliceJunctions.RDS")
