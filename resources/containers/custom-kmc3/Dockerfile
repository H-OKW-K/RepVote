FROM debian:bookworm-slim

ENV DEBIAN_FRONTEND=noninteractive

ARG TARGETARCH

RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    ca-certificates \
    gawk \
    coreutils \
    bash \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/kmc

RUN if [ "$TARGETARCH" = "arm64" ]; then \
        KMC_ARCH="arm64"; \
    else \
        KMC_ARCH="x64"; \
    fi && \
    echo "Downloading KMC3.2.4 for linux.${KMC_ARCH}..." && \
    mkdir -p /opt/kmc && \
    wget -q "https://github.com/refresh-bio/KMC/releases/download/v3.2.4/KMC3.2.4.linux.${KMC_ARCH}.tar.gz" -O /tmp/kmc.tar.gz && \
    tar -xzf /tmp/kmc.tar.gz -C /opt/kmc && \
    rm /tmp/kmc.tar.gz

ENV PATH="/opt/kmc/bin:/opt/kmc:${PATH}"

