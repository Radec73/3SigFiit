#FROM rocker/r-ver:4.3.2
FROM bioconductor/bioconductor_docker:RELEASE_3_18
MAINTAINER Radovan Cyprich

ENV R_VERSION=4.3.2
ENV PYTHON_VERSION=3.9.7
#ENV RENV_CONFIG_INSTALL_FROM="binary"
#ENV RENV_CONFIG_REPOS_OVERRIDE="https://cloud.r-project.org"
WORKDIR /app

RUN apt-get update && apt-get install -y software-properties-common && \
    apt-get update && apt-get install -y \
      apt-transport-https \
      python3-pip \
      python3 \
      libcurl4-openssl-dev \
      libssl-dev \
      libxml2-dev \
      liblzma-dev \
      libbz2-dev \
      libpcre2-dev \
      cmake \
      build-essential \
      curl \
      unzip \
      wget \
      sudo
      #r-base
RUN apt-get clean && rm -rf /var/lib/apt/lists/*

COPY . .

RUN R -e "install.packages('BiocManager', repos = 'https://cloud.r-project.org', ask = FALSE)" && \
    R -e "BiocManager::install(version = '3.18', ask = FALSE)" && \
    R -e "BiocManager::install(c('BSgenome', 'VariantAnnotation', 'BSgenome.Hsapiens.UCSC.hg19', 'BSgenome.Hsapiens.UCSC.hg38', 'maftools', 'sigminer', 'remotes', 'NMF', 'readr', 'IRanges', 'GenomicRanges', 'GenomeInfoDb', 'mutSignatures'), ask = FALSE, dependencies = TRUE)"


RUN pip install --upgrade pip setuptools wheel &&\
    pip install --upgrade pip --no-cache-dir -r requirements.txt

RUN python3 install_genome.py

CMD ["python3", "./superscript.py"]


