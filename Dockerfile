FROM rocker/r-ver:4.3.3

USER root

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    ca-certificates \
    curl \
    default-jre-headless \
    gfortran \
    git \
    libbz2-dev \
    libcurl4-openssl-dev \
    libgit2-dev \
    libglpk-dev \
    liblzma-dev \
    libncurses5-dev \
    libncursesw5-dev \
    libreadline-dev \
    libssl-dev \
    libxml2-dev \
    make \
    pkg-config \
    samtools \
    unzip \
    wget \
    xz-utils \
    zlib1g-dev && \
    rm -rf /var/lib/apt/lists/*

RUN mkdir -p /tmp/downloads

RUN curl -fsSL --retry 10 -o /tmp/downloads/htslib.tar.gz https://github.com/samtools/htslib/archive/1.7.tar.gz && \
    mkdir -p /tmp/downloads/htslib && \
    tar -C /tmp/downloads/htslib --strip-components 1 -zxf /tmp/downloads/htslib.tar.gz && \
    make -C /tmp/downloads/htslib && \
    rm -f /tmp/downloads/htslib.tar.gz

ENV HTSLIB=/tmp/downloads/htslib

RUN curl -fsSL --retry 10 -o /tmp/downloads/alleleCount.tar.gz https://github.com/cancerit/alleleCount/archive/v4.0.0.tar.gz && \
    mkdir -p /tmp/downloads/alleleCount && \
    tar -C /tmp/downloads/alleleCount --strip-components 1 -zxf /tmp/downloads/alleleCount.tar.gz && \
    mkdir -p /tmp/downloads/alleleCount/c/bin && \
    make -C /tmp/downloads/alleleCount/c && \
    cp /tmp/downloads/alleleCount/c/bin/alleleCounter /usr/local/bin/alleleCounter && \
    rm -rf /tmp/downloads/alleleCount /tmp/downloads/alleleCount.tar.gz

RUN curl -fsSL --retry 10 -o /tmp/downloads/impute2.tgz https://mathgen.stats.ox.ac.uk/impute/impute_v2.3.2_x86_64_static.tgz && \
    mkdir -p /tmp/downloads/impute2 && \
    tar -C /tmp/downloads/impute2 --strip-components 1 -zxf /tmp/downloads/impute2.tgz && \
    cp /tmp/downloads/impute2/impute2 /usr/local/bin/impute2 && \
    chmod +x /usr/local/bin/impute2 && \
    rm -rf /tmp/downloads/impute2 /tmp/downloads/impute2.tgz

RUN R -q -e 'install.packages(c("BiocManager"), repos="https://cloud.r-project.org")'

RUN R -q -e 'install.packages(c("RColorBrewer", "argparse", "data.table", "doParallel", "dplyr", "foreach", "ggplot2", "gridExtra", "gtools", "optparse", "parallel", "readr", "R.utils", "stringr", "tidyr"), repos="https://cloud.r-project.org")'

RUN mkdir -p /usr/local/lib/R/site-library && \
  R -q -e 'options(repos="https://cloud.r-project.org"); \
    install.packages("tidyverse", lib="/usr/local/lib/R/site-library", dependencies=TRUE)'
ENV R_LIBS_SITE=/usr/local/lib/R/site-library

RUN R -q -e 'BiocManager::install(c("GenomicRanges", "StructuralVariantAnnotation", "VariantAnnotation"), ask=FALSE, update=FALSE)'

RUN curl -fsSL --retry 10 -o /tmp/downloads/copynumber.tar.gz https://github.com/igordot/copynumber/archive/refs/heads/master.tar.gz && \
    R CMD INSTALL /tmp/downloads/copynumber.tar.gz && \
    rm -f /tmp/downloads/copynumber.tar.gz

RUN curl -fsSL --retry 10 -o /tmp/downloads/ascat.tar.gz https://github.com/VanLoo-lab/ascat/archive/refs/heads/master.tar.gz && \
    mkdir -p /tmp/downloads/ascat && \
    tar -C /tmp/downloads/ascat --strip-components 1 -zxf /tmp/downloads/ascat.tar.gz && \
    R CMD INSTALL /tmp/downloads/ascat/ASCAT && \
    rm -rf /tmp/downloads/ascat /tmp/downloads/ascat.tar.gz

ARG GITHUB_TOKEN=""
RUN if [ -z "$GITHUB_TOKEN" ]; then \
      git clone --depth 1 https://github.com/McGranahanLab/battenberg_wgs.git /opt/battenberg; \
    else \
      git clone --depth 1 https://${GITHUB_TOKEN}@github.com/McGranahanLab/battenberg_wgs.git /opt/battenberg; \
    fi

RUN R CMD INSTALL /opt/battenberg

CMD ["/bin/bash"]
