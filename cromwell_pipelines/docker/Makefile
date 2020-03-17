PWD=`pwd`

# tea

.PHONY : base-image base-scripts base-bam-tools base-vcf-tools base-other-tools str-callers-shared-tasks gangstr expansion-hunter expansion-hunter3 gatk stretch download-pacbio-reads repeat-finder simulate-locus evaluate-indel-calls

base-image:
	printf "Base Image \n" \
	&& cd $(PWD)/base-image && make | grep sha256 | grep latest | cut -d ' ' -f 3 > sha256.txt \
	&& cat sha256.txt && [ ! -z "`cat sha256.txt`" ] && sed -i.bak "s/base-image@sha256:[^\"]*/"base-image@`cat sha256.txt`"/"  ../../wdl/*/*.wdl \
	&& rm ../../wdl/*/*.wdl.bak \
	&& printf "============\n"

base-scripts:
	printf "Base Scripts Image \n" \
	&& cd $(PWD)/base-scripts && make | grep sha256 | grep latest | cut -d ' ' -f 3 > sha256.txt \
	&& printf "============\n"

base-bam-tools:
	printf "Base Bam Tools Image \n" \
	&& cd $(PWD)/base-bam-tools && make | grep sha256 | grep latest | cut -d ' ' -f 3 > sha256.txt \
	&& printf "============\n"

base-vcf-tools:
	printf "Base VCF Tools Image \n" \
	&& cd $(PWD)/base-vcf-tools && make | grep sha256 | grep latest | cut -d ' ' -f 3 > sha256.txt \
	&& printf "============\n"

base-other-tools:
	printf "Base Other Tools Image \n" \
	&& cd $(PWD)/base-other-tools && make | grep sha256 | grep latest | cut -d ' ' -f 3 > sha256.txt \
	&& printf "============\n"

gangstr: base-image
	printf "GangSTR \n" \
	&& cd $(PWD)/gangstr && make | grep sha256 | grep latest | cut -d ' ' -f 3 > sha256.txt \
	&& cat sha256.txt && [ ! -z "`cat sha256.txt`" ] && sed -i.bak "s/gangstr@sha256:[^\"]*/"gangstr@`cat sha256.txt`"/"  ../../wdl/*/*.wdl \
	&& rm ../../wdl/*/*.wdl.bak \
	&& printf "============\n"

expansion-hunter: base-image
	printf "ExpansionHunter \n" \
	&& cd $(PWD)/expansion-hunter && make | grep sha256 | grep latest | cut -d ' ' -f 3 > sha256.txt \
	&& cat sha256.txt && [ ! -z "`cat sha256.txt`" ] && sed -i.bak "s/expansion-hunter@sha256:[^\"]*/"expansion-hunter@`cat sha256.txt`"/"  ../../wdl/*/*.wdl \
	&& rm ../../wdl/*/*.wdl.bak \
	&& printf "============\n"

expansion-hunter3: base-image
	printf "ExpansionHunter3 \n" \
	&& cd $(PWD)/expansion-hunter3 && make | grep sha256 | grep latest | cut -d ' ' -f 3 > sha256.txt \
	&& cat sha256.txt && [ ! -z "`cat sha256.txt`" ] && sed -i.bak "s/expansion-hunter3@sha256:[^\"]*/"expansion-hunter3@`cat sha256.txt`"/"  ../../wdl/*/*.wdl \
	&& rm ../../wdl/*/*.wdl.bak \
	&& printf "============\n"

gatk: base-image
	printf "Gatk \n" \
	&& cd $(PWD)/gatk && make | grep sha256 | grep latest | cut -d ' ' -f 3 > sha256.txt \
	&& cat sha256.txt && [ ! -z "`cat sha256.txt`" ] && sed -i.bak "s/gatk@sha256:[^\"]*/"gatk@`cat sha256.txt`"/"  ../../wdl/*/*.wdl \
	&& rm ../../wdl/*/*.wdl.bak \
	&& printf "============\n"

stretch: base-image
	printf "Stretch \n" \
	&& cd $(PWD)/stretch && make | grep sha256 | grep latest | cut -d ' ' -f 3 > sha256.txt \
	&& cat sha256.txt && [ ! -z "`cat sha256.txt`" ] && sed -i.bak "s/stretch@sha256:[^\"]*/"stretch@`cat sha256.txt`"/"  ../../wdl/*/*.wdl \
	&& rm ../../wdl/*/*.wdl.bak \
	&& printf "============\n"

tredparse: base-image
	printf "Tredparse \n" \
	&& cd $(PWD)/tredparse && make | grep sha256 | grep latest | cut -d ' ' -f 3 > sha256.txt \
	&& cat sha256.txt && [ ! -z "`cat sha256.txt`" ] && sed -i.bak "s/tredparse@sha256:[^\"]*/"tredparse@`cat sha256.txt`"/"  ../../wdl/*/*.wdl \
	&& rm ../../wdl/*/*.wdl.bak \
	&& printf "============\n"

download-pacbio-reads: base-image
	printf "Download PacBio Reads \n" \
	&& cd $(PWD)/download-pacbio-reads && make | grep sha256 | grep latest | cut -d ' ' -f 3 > sha256.txt \
	&& cat sha256.txt && [ ! -z "`cat sha256.txt`" ] && sed -i.bak "s/download-pacbio-reads@sha256:[^\"]*/"download-pacbio-reads@`cat sha256.txt`"/"  ../../wdl/*/*.wdl \
	&& rm ../../wdl/*/*.wdl.bak \
	&& printf "============\n"

repeat-finder: base-image
	printf "Repeat Finder \n" \
	&& cd $(PWD)/repeat-finder && make | grep sha256 | grep latest | cut -d ' ' -f 3 > sha256.txt \
	&& cat sha256.txt && [ ! -z "`cat sha256.txt`" ] && sed -i.bak "s/repeat-finder@sha256:[^\"]*/"repeat-finder@`cat sha256.txt`"/"  ../../wdl/*/*.wdl \
	&& rm ../../wdl/*/*.wdl.bak \
	&& printf "============\n"

tea: base-image
	printf "TEA \n" \
	&& cd $(PWD)/tea && make | grep sha256 | grep latest | cut -d ' ' -f 3 > sha256.txt \
	&& cat sha256.txt && [ ! -z "`cat sha256.txt`" ] && sed -i.bak "s/tea@sha256:[^\"]*/"tea@`cat sha256.txt`"/"  ../../wdl/*/*.wdl \
	&& rm ../../wdl/*/*.wdl.bak \
	&& printf "============\n"

simulate-locus: base-image base-scripts base-bam-tools
	printf "Simulate Locus \n" \
	&& cd $(PWD)/simulate-locus && make | grep sha256 | grep latest | cut -d ' ' -f 3 > sha256.txt \
	&& cat sha256.txt && [ ! -z "`cat sha256.txt`" ] && sed -i.bak "s/simulate-locus@sha256:[^\"]*/"simulate-locus@`cat sha256.txt`"/"  ../../wdl/*/*.wdl \
	&& rm ../../wdl/*/*.wdl.bak \
	&& printf "============\n"

evaluate-indel-calls: base-image base-scripts base-bam-tools
	printf "Evaluate InDel Calls \n" \
	&& cd $(PWD)/evaluate-indel-calls && make | grep sha256 | grep latest | cut -d ' ' -f 3 > sha256.txt \
	&& cat sha256.txt && [ ! -z "`cat sha256.txt`" ] && sed -i.bak "s/evaluate-indel-calls@sha256:[^\"]*/"evaluate-indel-calls@`cat sha256.txt`"/"  ../../wdl/*/*.wdl \
	&& rm ../../wdl/*/*.wdl.bak \
	&& printf "============\n"

str-callers-shared-tasks: base-image base-scripts base-bam-tools base-vcf-tools base-other-tools
	printf "STR Callers' Shared Tasks \n" \
	&& cd $(PWD)/str-callers-shared-tasks && make | grep sha256 | grep latest | cut -d ' ' -f 3 > sha256.txt \
	&& cat sha256.txt && [ ! -z "`cat sha256.txt`" ] && sed -i.bak "s/str-callers-shared-tasks@sha256:[^\"]*/"str-callers-shared-tasks@`cat sha256.txt`"/"  ../../wdl/*/*.wdl \
	&& rm ../../wdl/*/*.wdl.bak \
	&& printf "============\n"

all: base-image base-scripts base-bam-tools base-vcf-tools base-other-tools str-callers-shared-tasks evaluate-indel-calls simulate-locus gangstr expansion-hunter expansion-hunter3 gatk stretch
	#tredparse
	#repeat-finder
	#docker system prune -f
