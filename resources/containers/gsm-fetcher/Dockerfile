# Use a lightweight, standard Linux base
FROM ubuntu:22.04

# Prevent interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# 1. Install system dependencies
# - wget/curl: for downloading NCBI tools
# - pigz: for the 'dump_merge_and_layout' rule (parallel gzip)
# - uuid-runtime: required for generating SRA config
# - python3: often required by internal SRA scripts or Snakemake wrappers
# - perl: required by Entrez Direct
RUN apt-get update && apt-get install -y \
    wget \
    curl \
    pigz \
    python3 \
    perl \
    uuid-runtime \
    libxml2 \
    && rm -rf /var/lib/apt/lists/*

# 2. Install SRA Toolkit (prefetch, fasterq-dump)
# We use the official pre-compiled Ubuntu binary to ensure stability
ENV SRA_VERSION=3.1.1
WORKDIR /opt
RUN wget "https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/${SRA_VERSION}/sratoolkit.${SRA_VERSION}-ubuntu64.tar.gz" \
    && tar -xzf sratoolkit.${SRA_VERSION}-ubuntu64.tar.gz \
    && rm sratoolkit.${SRA_VERSION}-ubuntu64.tar.gz

# Add SRA Toolkit to PATH
ENV PATH="/opt/sratoolkit.${SRA_VERSION}-ubuntu64/bin:${PATH}"

# Configure SRA Toolkit
# This step is critical to prevent 'fasterq-dump' from failing in a container environment
# It generates a basic configuration file
RUN mkdir -p /root/.ncbi && \
    printf '/LIBS/GUID = "%s"\n' `uuidgen` > /root/.ncbi/user-settings.mkfg && \
    vdb-config --restore-defaults

# 3. Install Entrez Direct (esearch, efetch)
# Installs into /root/edirect
WORKDIR /root
RUN sh -c "$(curl -fsSL https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"

# Add Entrez Direct to PATH
ENV PATH="/root/edirect:${PATH}"

# 4. Final Setup
WORKDIR /data
CMD ["/bin/bash"]
#docker build -t gsm-fetcher:1.0.0 . && apptainer build gsm-fetcher.sif docker-daemon://gsm-fetcher:1.0.0
