FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    wget \
    curl \
    pigz \
    python3 \
    perl \
    uuid-runtime \
    libxml2 \
    && rm -rf /var/lib/apt/lists/*

ENV SRA_VERSION=3.1.1
WORKDIR /opt
RUN wget "https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/${SRA_VERSION}/sratoolkit.${SRA_VERSION}-ubuntu64.tar.gz" \
    && tar -xzf sratoolkit.${SRA_VERSION}-ubuntu64.tar.gz \
    && rm sratoolkit.${SRA_VERSION}-ubuntu64.tar.gz

ENV PATH="/opt/sratoolkit.${SRA_VERSION}-ubuntu64/bin:${PATH}"

RUN mkdir -p /root/.ncbi && \
    printf '/LIBS/GUID = "%s"\n' `uuidgen` > /root/.ncbi/user-settings.mkfg && \
    vdb-config --restore-defaults

WORKDIR /root
RUN sh -c "$(curl -fsSL https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"

ENV PATH="/root/edirect:${PATH}"

WORKDIR /data
CMD ["/bin/bash"]

