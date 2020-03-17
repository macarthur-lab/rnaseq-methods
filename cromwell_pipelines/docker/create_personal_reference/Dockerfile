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
    liblzma-dev

# install RScript
RUN install_packages libopenblas-base  r-base

# install java-8    (adds 340Mb)
RUN install_packages default-jdk

# install python2 utils
RUN install_packages python-dev python-setuptools python-pip \
	&& pip install --upgrade wheel

# install python3.7 (adds 800Mb)
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
RUN python3 -m pip install --upgrade pip \
	&& python3 -m pip install --upgrade tqdm pyvcf pyfaidx pybedtools cyvcf2 peewee scipy numpy pandas matplotlib seaborn ipython

# nice-to-have linux utils
RUN install_packages less htop vim emacs

COPY bashrc /root/.bashrc


