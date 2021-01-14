FROM bitnami/minideb:stretch

MAINTAINER Ben Weisburd

RUN install_packages \
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
    liblzma-dev \
    python3-dev \
    git

# install gcloud
ENV GCLOUD_SDK_VERISON="321.0.0"
RUN wget https://dl.google.com/dl/cloudsdk/channels/rapid/downloads/google-cloud-sdk-${GCLOUD_SDK_VERISON}-linux-x86_64.tar.gz \
    && tar xzf google-cloud-sdk-${GCLOUD_SDK_VERISON}-linux-x86_64.tar.gz \
    && /google-cloud-sdk/install.sh -q \
    && /google-cloud-sdk/bin/gcloud components update --quiet

ENV PATH=/google-cloud-sdk/bin:$PATH

RUN wget https://github.com/brentp/somalier/files/4566475/sites.hg38.rna.vcf.gz
RUN wget https://github.com/brentp/somalier/releases/download/v0.2.12/somalier -O /usr/bin/somalier \
    && chmod 777 /usr/bin/somalier

ENTRYPOINT /bin/bash
