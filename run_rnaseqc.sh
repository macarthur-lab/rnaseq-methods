#!/bin/bash

BAM=$1
GTF_REFERENCE=$2
SAMPLE_NAME=$(basename ${BAM})

# run rnaseqc
	rnaseqc.v2.1.2.linux ${GTF_REFERENCE} ${BAM} --coverage .
# 

if [ ! -f ${BAM}.gene_tpm.gct ]
	exit
fi

tail -n+4 ${BAM}.gene_tpm.gct | sort | cut -f3 >> gene_tpm.gct

