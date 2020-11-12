FROM bitnami/minideb:buster

MAINTAINER Ben Weisburd

RUN install_packages \
    ca-certificates \
    apt-utils \
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
    libssl-dev \
    libcurl4-gnutls-dev \
    libbz2-dev \
    liblzma-dev \
    libxml2-dev \
    gfortran \
    libreadline-dev \
    libpcre3 \
    libpcre3-dev \
    libpng-dev \
    libjpeg-dev \
    default-libmysqld-dev \
    git \
    libcairo2-dev

#RUN install_packages default-jdk

# install RScript

RUN wget https://cran.r-project.org/src/base/R-4/R-4.0.0.tar.gz \
    && tar xzf R-4.0.0.tar.gz \
    && rm R-4.0.0.tar.gz \
    && cd R-4.0.0 \
    && ./configure --with-x=no --with-pcre1 --with-cairo=yes --with-libpng=yes \
    && make \
    && make install

## install python3.7 (adds 800Mb)
RUN apt-get update \
	&& apt-get dist-upgrade -y \
	&& apt-get install -y python-smbus libncursesw5-dev libgdbm-dev libc6-dev zlib1g-dev libsqlite3-dev libssl-dev openssl libffi-dev

ENV PYTHON_VERSION="3.7.1"
RUN wget https://www.python.org/ftp/python/${PYTHON_VERSION}/Python-${PYTHON_VERSION}.tar.xz \
	&& tar xvf Python-${PYTHON_VERSION}.tar.xz \
	&& rm Python-${PYTHON_VERSION}.tar.xz \
	&& cd Python-${PYTHON_VERSION} \
	&& ./configure --enable-optimizations \
	&& make install \
	&& make clean

# python3 packages
RUN python3 -m pip install --upgrade pip

# nice-to-have linux utils
RUN install_packages less htop vim emacs procps

# install R pacakges
RUN R -e "install.packages('BiocManager', repos='http://cran.rstudio.com')"
#RUN R -e "BiocManager::install('FRASER')"
RUN R -e "install.packages('devtools', repos='http://cran.rstudio.com')"
RUN R -e "devtools::install_github('gagneurlab/FRASER', dependencies=TRUE)"
RUN R -e "install.packages('argparse', repos='http://cran.rstudio.com')"
RUN R -e "BiocManager::install('OUTRIDER')"
RUN R -e "BiocManager::install('mumichae/tMAE')"
RUN R -e "BiocManager::install('stephenturner/annotables')"
RUN R -e "BiocManager::install('Biostrings')"
RUN R -e "BiocManager::install('IRanges')"
RUN R -e "BiocManager::install('DelayedArray')"
RUN R -e "install.packages('data.table', repos='http://cran.rstudio.com')"
RUN R -e "install.packages('ggplot2', repos='http://cran.rstudio.com')"
RUN R -e "install.packages('ggpubr', repos='http://cran.rstudio.com')"
RUN R -e "install.packages('dplyr', repos='http://cran.rstudio.com')"
RUN R -e "install.packages('purrr', repos='http://cran.rstudio.com')"
RUN R -e "install.packages('ggrepel', repos='http://cran.rstudio.com')"
RUN R -e "install.packages('plotly', repos='http://cran.rstudio.com')"
RUN R -e "install.packages('stringr', repos='http://cran.rstudio.com')"

RUN install_packages libxt-dev
RUN R -e "install.packages('Cairo', repos='http://cran.rstudio.com')"
RUN install_packages xvfb xauth

# install gcloud
ENV GCLOUD_SDK_VERISON="299.0.0"
RUN wget https://dl.google.com/dl/cloudsdk/channels/rapid/downloads/google-cloud-sdk-${GCLOUD_SDK_VERISON}-linux-x86_64.tar.gz \
    && tar xzf google-cloud-sdk-${GCLOUD_SDK_VERISON}-linux-x86_64.tar.gz \
    && /google-cloud-sdk/install.sh -q \
    && /google-cloud-sdk/bin/gcloud components update --quiet \
    && rm -rf /google-cloud-sdk/.install/

#RUN R -e "devtools::install_github('gagneurlab/FRASER', ref='R3.6', dependencies=TRUE)"
#RUN R -e "devtools::install_github('gagneurlab/OUTRIDER', dependencies=TRUE)"

RUN install_packages bc graphviz netbase gnupg

# install java-8    (adds 340Mb)
RUN install_packages default-jdk

# install gatk
ENV GATK_VERSION="4.1.9.0"
RUN wget https://github.com/broadinstitute/gatk/releases/download/${GATK_VERSION}/gatk-${GATK_VERSION}.zip \
    && unzip gatk-${GATK_VERSION}.zip

ENV PATH=/gatk-${GATK_VERSION}:$PATH

# install htslib
ENV SAMTOOLS_VERSION="1.11"
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

# install samtools
RUN wget https://github.com/samtools/bcftools/releases/download/${SAMTOOLS_VERSION}/bcftools-${SAMTOOLS_VERSION}.tar.bz2 \
	&& tar xjf bcftools-${SAMTOOLS_VERSION}.tar.bz2 \
	&& rm bcftools-${SAMTOOLS_VERSION}.tar.bz2 \
	&& cd bcftools-${SAMTOOLS_VERSION} \
	&& ./configure \
	&& make \
	&& make install \
	&& make clean


# install pandoc
ENV PANDOC_VERSION="1.17.0.3"
RUN wget https://hackage.haskell.org/package/pandoc-${PANDOC_VERSION}/pandoc-${PANDOC_VERSION}.tar.gz \
    && tar xvzf pandoc-${PANDOC_VERSION}.tar.gz \
    && rm pandoc-${PANDOC_VERSION}.tar.gz \
    && cd pandoc-${PANDOC_VERSION} \
    && curl -sSL https://get.haskellstack.org/ | sh \
    && stack setup \
    && stack install \
    && rm -rf /root/.stack/

ENV PATH=/root/.local/bin:$PATH

# install drop
RUN python3 -m pip install --upgrade pip
RUN python3 -m pip install git+https://github.com/gagneurlab/drop.git@1.0.1
RUN python3 -m pip install snakemake
RUN mkdir ~/drop_demo \
    && cd ~/drop_demo \
    && drop demo \
    && snakemake --list-target-rules \
    && rm -rf resource/rna_bam/ resource/dna_vcf/ resource/*.vcf.gz resource/chr21.fa.gz \
    && rm /tmp/resource.tar.gz

ENV PATH=/google-cloud-sdk/bin:$PATH

#COPY *.R /
#COPY .Rprofile /root/
#COPY bashrc /root/.bashrc

WORKDIR /
