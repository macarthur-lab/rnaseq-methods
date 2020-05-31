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
    git

#RUN install_packages default-jdk

# install RScript

RUN wget https://cran.r-project.org/src/base/R-4/R-4.0.0.tar.gz \
    && tar xzf R-4.0.0.tar.gz \
    && rm R-4.0.0.tar.gz \
    && cd R-4.0.0 \
    && ./configure --with-x=no --with-pcre1 \
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

# install gcloud
ENV GCLOUD_SDK_VERISON="293.0.0"
RUN wget https://dl.google.com/dl/cloudsdk/channels/rapid/downloads/google-cloud-sdk-${GCLOUD_SDK_VERISON}-linux-x86_64.tar.gz \
    && tar xzf google-cloud-sdk-${GCLOUD_SDK_VERISON}-linux-x86_64.tar.gz \
    && /google-cloud-sdk/install.sh -q \
    && /google-cloud-sdk/bin/gcloud components update --quiet


# nice-to-have linux utils
RUN install_packages less htop vim emacs procps

# install R pacakges
RUN R -e "install.packages('BiocManager', repos='http://cran.rstudio.com')"
#RUN R -e "BiocManager::install('FRASER')"
RUN R -e "install.packages('devtools', repos='http://cran.rstudio.com')"
RUN R -e "devtools::install_github('gagneurlab/FRASER', dependencies=TRUE)"
RUN R -e "install.packages('argparse', repos='http://cran.rstudio.com')"

#RUN R -e "BiocManager::install('OUTRIDER')"
#RUN R -e "BiocManager::install('mumichae/tMAE')"

#RUN R -e "devtools::install_github('gagneurlab/FRASER', ref='R3.6', dependencies=TRUE)"
#RUN R -e "devtools::install_github('gagneurlab/OUTRIDER', dependencies=TRUE)"

ENV PATH=/google-cloud-sdk/bin:$PATH

COPY countRNA.R /

WORKDIR /

#COPY bashrc /root/.bashrc
