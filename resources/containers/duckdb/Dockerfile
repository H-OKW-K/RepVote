FROM python:3.11-slim-bookworm

RUN apt-get update && apt-get install -y --no-install-recommends \
    gawk \
    wget \
    unzip \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

RUN wget https://github.com/duckdb/duckdb/releases/download/v1.1.3/duckdb_cli-linux-amd64.zip \
    && unzip duckdb_cli-linux-amd64.zip -d /usr/local/bin \
    && chmod +x /usr/local/bin/duckdb \
    && rm duckdb_cli-linux-amd64.zip

RUN pip install --no-cache-dir \
    duckdb==1.1.3 \
    pandas \
    pyarrow

#Verification
RUN gawk --version \
    && duckdb --version \
    && python3 -c "import duckdb; print('DuckDB Python lib works')"

CMD ["python3"]
