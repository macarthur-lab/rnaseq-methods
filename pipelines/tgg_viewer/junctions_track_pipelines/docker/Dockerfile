FROM python:3.8-slim

MAINTAINER Ben Weisburd


RUN apt-get update \
  && apt-get install -y \
    python-smbus \
    libncursesw5-dev \
    libgdbm-dev \
    libc6-dev \
    zlib1g-dev \
    libsqlite3-dev \
    libssl-dev \
    openssl \
    libffi-dev \
    ca-certificates \
    wget \
    curl \
    pkg-config \
    gcc \
    man-db \
    g++ \
    make \
    autoconf \
    unzip \
    bzip2 \
    zlib1g-dev \
    dpkg-dev \
    build-essential \
    libcurl4-openssl-dev \
    libbz2-dev \
    liblzma-dev

# install gcloud
ENV GCLOUD_SDK_VERISON="293.0.0"
RUN wget https://dl.google.com/dl/cloudsdk/channels/rapid/downloads/google-cloud-sdk-${GCLOUD_SDK_VERISON}-linux-x86_64.tar.gz \
    && tar xzf google-cloud-sdk-${GCLOUD_SDK_VERISON}-linux-x86_64.tar.gz \
    && /google-cloud-sdk/install.sh -q \
    && /google-cloud-sdk/bin/gcloud components update --quiet

ENV PATH=/google-cloud-sdk/bin:$PATH

RUN wget -P /usr/local/bin http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigWigMerge \
    && chmod 777 /usr/local/bin/bigWigMerge

RUN wget -P /usr/local/bin http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig \
    && chmod 777 /usr/local/bin/bedGraphToBigWig


# install htslib  --  install in downstream images with:
#     COPY --from=weisburd/base-bam-tools:latest  /usr/local/bin/tabix /usr/local/bin/tabix
#     COPY --from=weisburd/base-bam-tools:latest  /usr/local/bin/htsfile /usr/local/bin/htsfile
#     COPY --from=weisburd/base-bam-tools:latest  /usr/local/bin/bgzip /usr/local/bin/bgzip
#     COPY --from=weisburd/base-bam-tools:latest  /usr/local/bin/samtools /usr/local/bin/samtools
ENV SAMTOOLS_VERSION="1.12"
RUN wget https://github.com/samtools/htslib/releases/download/${SAMTOOLS_VERSION}/htslib-${SAMTOOLS_VERSION}.tar.bz2 \
	&& tar xjf htslib-${SAMTOOLS_VERSION}.tar.bz2 \
	&& rm htslib-${SAMTOOLS_VERSION}.tar.bz2 \
	&& cd htslib-${SAMTOOLS_VERSION} \
	&& ./configure \
	&& make \
	&& make install \
	&& make clean


COPY gencode.v26.annotation.gff3.gz /
COPY gencode.v35.annotation.gff3.gz /

RUN python3 -m pip install pandas psutil

COPY combine_splice_junctions_using_pandas.py /
COPY convert_SJ_out_tab_to_junctions_bed.py /
COPY gencode_utils.py /

ENTRYPOINT /bin/bash

