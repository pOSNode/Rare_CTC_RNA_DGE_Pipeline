# ============================================================
# Dockerfile: ctc-rnaseq pipeline
# Single image containing all tools required for the pipeline.
# Base: rocker/r-ver for reproducible R environment.
# ============================================================

FROM rocker/r-ver:4.3.2

LABEL maintainer="Owais Siddiqi"
LABEL description="CTC low-input RNA-seq pipeline: fastp, Salmon, MultiQC, edgeR"
LABEL version="1.0.0"

# System dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    curl \
    unzip \
    python3 \
    python3-pip \
    bc \
    && rm -rf /var/lib/apt/lists/*

# ------------------------------------------------------------
# fastp (v0.23.4)
# ------------------------------------------------------------
RUN wget -q https://github.com/OpenGene/fastp/releases/download/v0.23.4/fastp \
    -O /usr/local/bin/fastp && \
    chmod +x /usr/local/bin/fastp

# ------------------------------------------------------------
# Salmon (v1.10.3)
# ------------------------------------------------------------
RUN wget -q https://github.com/COMBINE-lab/salmon/releases/download/v1.10.3/salmon-1.10.3_linux_x86_64.tar.gz \
    -O /tmp/salmon.tar.gz && \
    tar -xzf /tmp/salmon.tar.gz -C /opt && \
    ln -s /opt/salmon-1.10.3_linux_x86_64/bin/salmon /usr/local/bin/salmon && \
    rm /tmp/salmon.tar.gz

# ------------------------------------------------------------
# MultiQC
# ------------------------------------------------------------
RUN pip3 install --no-cache-dir multiqc==1.21

# ------------------------------------------------------------
# R packages: tximport, edgeR, limma, annotation
# ------------------------------------------------------------
RUN Rscript -e " \
    install.packages('BiocManager', repos='https://cloud.r-project.org'); \
    BiocManager::install(c( \
        'tximport', \
        'edgeR', \
        'limma', \
        'org.Hs.eg.db', \
        'AnnotationDbi', \
        'GenomicFeatures' \
    ), ask=FALSE); \
    install.packages(c( \
        'optparse', \
        'data.table', \
        'tidyverse' \
    ), repos='https://cloud.r-project.org')"

# ------------------------------------------------------------
# Copy pipeline R scripts into image
# ------------------------------------------------------------
COPY scripts/ /scripts/
RUN chmod +x /scripts/*.R

# Sanity checks at build time
RUN fastp --version 2>&1 | head -1
RUN salmon --version 2>&1 | head -1
RUN multiqc --version
RUN Rscript -e "library(edgeR); library(tximport); cat('R packages OK\n')"

WORKDIR /data
