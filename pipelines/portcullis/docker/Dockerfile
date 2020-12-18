FROM bitnami/minideb:buster

MAINTAINER Ben Weisburd

RUN apt update \
  && apt install -y \
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
ENV GCLOUD_SDK_VERISON="321.0.0"
RUN wget https://dl.google.com/dl/cloudsdk/channels/rapid/downloads/google-cloud-sdk-${GCLOUD_SDK_VERISON}-linux-x86_64.tar.gz \
    && tar xzf google-cloud-sdk-${GCLOUD_SDK_VERISON}-linux-x86_64.tar.gz \
    && /google-cloud-sdk/install.sh -q \
    && /google-cloud-sdk/bin/gcloud components update --quiet

ENV PATH=/google-cloud-sdk/bin:$PATH

# install htslib
ENV SAMTOOLS_VERSION="1.9"
RUN wget https://github.com/samtools/htslib/releases/download/${SAMTOOLS_VERSION}/htslib-${SAMTOOLS_VERSION}.tar.bz2 \
	&& tar xjf htslib-${SAMTOOLS_VERSION}.tar.bz2 \
	&& rm htslib-${SAMTOOLS_VERSION}.tar.bz2 \
	&& cd htslib-${SAMTOOLS_VERSION} \
	&& ./configure \
	&& make \
	&& make install \
	&& make clean

# install samtools
RUN wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 \
	&& tar xjf samtools-${SAMTOOLS_VERSION}.tar.bz2 \
	&& rm samtools-${SAMTOOLS_VERSION}.tar.bz2 \
	&& cd samtools-${SAMTOOLS_VERSION} \
	&& ./configure \
	&& make \
	&& make install \
	&& make clean

#RUN python3 -m pip install pandas
#RUN python -m pip install pandas
#ENV PYTHONPATH=/usr/local/lib/python3.7/dist-packages/

RUN apt update && apt-get install -y python-pip python-pandas python3-pip python3-pandas
RUN apt update && apt-get install -y git libtool libboost-all-dev

RUN git clone https://github.com/maplesond/portcullis.git \
    && cd portcullis \
    && ./autogen.sh \
    && ./configure \
    && make -j 4 \
    && make install


ENTRYPOINT /bin/bash

