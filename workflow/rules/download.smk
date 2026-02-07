# workflow/rules/download.smk

import os
import shutil
import pandas as pd

# ------------------------------------------------------------------------------
# 1. Helper Functions
# ------------------------------------------------------------------------------

def get_gsm_id(wildcards):
    """Retrieve the GSM ID from the samples dataframe."""
    return samples_df.loc[wildcards.sample, "r1"]

def get_local_paths(wildcards):
    """Retrieve raw file paths for local samples."""
    row = samples_df.loc[wildcards.sample]
    return {"r1": row["r1"], "r2": row["r2"] if pd.notna(row["r2"]) else None}

def get_local_input_files(wildcards):
    """
    Returns a list of local files to ensure Snakemake checks 
    their existence before running the rule.
    """
    row = samples_df.loc[wildcards.sample]
    files = [row["r1"]]
    if pd.notna(row["r2"]) and str(row["r2"]).strip():
        files.append(row["r2"])
    return files

def get_sample_type(wildcards):
    """Safely retrieves the sample type (gsm or local)."""
    return samples_df.loc[wildcards.sample, "type"].strip().lower()

def get_staging_r1(wildcards):
    """Returns download rule output for GSM, or filepath for Local."""
    stype = get_sample_type(wildcards)
    if stype == "gsm":
        # Dependency: The output of the GSM processing rule
        return f"data/{wildcards.sample}/raw_gsm/reads_1.fastq.gz"
    else:
        # Dependency: Local file path on disk
        return samples_df.loc[wildcards.sample, "r1"]

def get_staging_r2(wildcards):
    """Returns download rule output for GSM, or filepath for Local (if exists)."""
    stype = get_sample_type(wildcards)
    if stype == "gsm":
        return f"data/{wildcards.sample}/raw_gsm/reads_2.fastq.gz"
    else:
        r2 = samples_df.loc[wildcards.sample, "r2"]
        # Return list if present, else empty list (Snakemake standard for optional inputs)
        return [r2] if (pd.notna(r2) and str(r2).strip()) else []

def get_staging_layout(wildcards):
    """Returns layout file for GSM (Local layout is generated on the fly)."""
    if get_sample_type(wildcards) == "gsm":
        return [f"data/{wildcards.sample}/raw_gsm/layout.txt"]
    return []

# ------------------------------------------------------------------------------
# 2. GSM Download Pipeline
# ------------------------------------------------------------------------------

rule find_srrs:
    """Find associated SRR IDs for a given GSM ID using esearch."""
    container: config["docker_images"]["gsm_fetcher"]
    output:
        srr_list = os.path.join(config["paths"]["data"], "{sample}/raw_gsm/srr_ids.txt")
    params:
        gsm_id = get_gsm_id
    log:
        os.path.join(config["paths"]["data"], "{sample}/logs/find_srrs.log")
    resources:
        io_bound = config["resource_profiles"]["io_bound_for_esearch"]
    shell:
        """
        echo "Finding SRRs for {params.gsm_id}..." > "{log}"
        
        # Use NCBI E-utilities to fetch run info
        esearch -db sra -query "{params.gsm_id}" | \
        efetch -format runinfo | \
        cut -d',' -f1 | \
        grep "SRR" > {output.srr_list} 2>> "{log}"

        if [ ! -s {output.srr_list} ]; then
            echo "Error: No SRR files found for {params.gsm_id}." >> "{log}"
            exit 1
        fi
        """

rule prefetch_sra:
    """Download .sra files using prefetch."""
    container: config["docker_images"]["gsm_fetcher"]
    input:
        srr_list = os.path.join(config["paths"]["data"], "{sample}/raw_gsm/srr_ids.txt")
    output:
        sra_dir = temp(directory(os.path.join(config["paths"]["data"], "{sample}/raw_gsm/sra_files")))
    log:
        os.path.join(config["paths"]["data"], "{sample}/logs/prefetch.log")
    resources: 
        io_bound = config["resource_profiles"]["io_bound_for_prefetch"]
    shell:
        """
        mkdir -p {output.sra_dir}
        echo "Prefetching SRR files..." > "{log}"
        
        while read srr_id || [ -n "$srr_id" ]; do
            if [ -n "$srr_id" ]; then
                echo "Prefetching $srr_id..." >> "{log}"
                # --max-size 100G allows large files
                prefetch "$srr_id" -O {output.sra_dir} --max-size 100G >> "{log}" 2>&1
            fi
        done < {input.srr_list}
        """

rule dump_merge_and_layout_gsm:
    """
    Dump SRA to Fastq, merge multiple runs if necessary, and detect layout (PE/SE).
    Outputs standard named files (reads_1, reads_2) to simplify downstream rules.
    """
    container: config["docker_images"]["gsm_fetcher"]
    input:
        srr_list = os.path.join(config["paths"]["data"], "{sample}/raw_gsm/srr_ids.txt"),
        sra_dir = os.path.join(config["paths"]["data"], "{sample}/raw_gsm/sra_files")
    output:
        r1 = os.path.join(config["paths"]["data"], "{sample}/raw_gsm/reads_1.fastq.gz"),
        r2 = os.path.join(config["paths"]["data"], "{sample}/raw_gsm/reads_2.fastq.gz"),
        layout = os.path.join(config["paths"]["data"], "{sample}/raw_gsm/layout.txt")
    params:
        temp_dir = os.path.join(config["paths"]["data"], "{sample}/raw_gsm/temp_fastq")
    log:
        os.path.join(config["paths"]["data"], "{sample}/logs/dump_merge.log")
    threads: config["resource_profiles"]["heavy"]["threads"]
    resources:
        mem_mb = config["resource_profiles"]["heavy"]["mem_mb"]
    shell:
        """
        # 0. Setup
        mkdir -p {params.temp_dir}
        echo "Starting fastq-dump..." > "{log}"

        # 1. Dump all SRRs associated with this GSM
        while read srr_id || [ -n "$srr_id" ]; do
            if [ -n "$srr_id" ]; then
                echo "Dumping $srr_id..." >> "{log}"
                fasterq-dump "{input.sra_dir}/$srr_id" \
                    --split-files \
                    --outdir {params.temp_dir} \
                    --threads {threads} \
                    >> "{log}" 2>&1
            fi
        done < {input.srr_list}

        # 2. Detect Layout & Merge
        # Check if we have any *_2.fastq files indicating Paired-End
        if ls {params.temp_dir}/*_2.fastq >/dev/null 2>&1; then
            echo "paired" > {output.layout}
            echo "Paired-end detected. Merging..." >> "{log}"
            
            # Merge all R1s and compress
            cat {params.temp_dir}/*_1.fastq | pigz -p {threads} > {output.r1}
            # Merge all R2s and compress
            cat {params.temp_dir}/*_2.fastq | pigz -p {threads} > {output.r2}
        else
            echo "single" > {output.layout}
            echo "Single-end detected. Merging..." >> "{log}"
            
            # Merge all files (should be just .fastq) into R1
            cat {params.temp_dir}/*.fastq | pigz -p {threads} > {output.r1}
            
            # Create a dummy R2 file to satisfy Snakemake's output requirement
            touch {output.r2}
        fi

        # 3. Cleanup
        rm -rf {params.temp_dir}
        echo "Dump and merge complete." >> "{log}"
        """

# ------------------------------------------------------------------------------
# 3. Unified File Staging
# ------------------------------------------------------------------------------

rule stage_combined_data:
    """
    Unifies Local and GSM data into a single 'raw' directory structure.
    - If GSM: Symlinks the downloaded/merged files from raw_gsm.
    - If Local: Symlinks the provided paths and generates a layout.txt.
    """
    container: config["docker_images"]["base_unix"]
    input:
        r1 = get_staging_r1,
        r2 = get_staging_r2,
        layout = get_staging_layout
    output:
        r1 = os.path.join(config["paths"]["data"], "{sample}/raw/reads_1.fastq.gz"),
        r2 = os.path.join(config["paths"]["data"], "{sample}/raw/reads_2.fastq.gz"),
        layout = os.path.join(config["paths"]["data"], "{sample}/raw/layout.txt")
    params:
        stype = get_sample_type
    log:
        os.path.join(config["paths"]["data"], "{sample}/logs/stage_combined.log")
    shell:
        """
        echo "Staging data for type: {params.stype}" > "{log}"

        # Helper function to safely link or copy
        safe_stage() {{
            local src_raw="$1"
            local dest="$2"
            
            # Resolve absolute paths to check if they are the same
            local src=$(readlink -f "$src_raw")
            local target=$(readlink -f "$dest" 2>/dev/null || echo "$dest")

            if [ "$src" != "$target" ]; then
                echo "Staging $src to $dest" >> "{log}"
                # Try symlink; if it fails (Operation not permitted), copy the file
                ln -sf "$src" "$dest" || cp "$src" "$dest"
            else
                echo "Source and destination are the same ($src), skipping." >> "{log}"
            fi
        }}

        if [ "{params.stype}" == "gsm" ]; then
            # --- GSM STRATEGY ---
            safe_stage "{input.r1}" "{output.r1}"
            safe_stage "{input.r2}" "{output.r2}"
            cp {input.layout} {output.layout}
            echo "Staged GSM files." >> "{log}"

        else
            # --- LOCAL STRATEGY ---
            safe_stage "{input.r1}" "{output.r1}"

            if [ -n "{input.r2}" ]; then
                safe_stage "{input.r2}" "{output.r2}"
                echo "paired" > {output.layout}
                echo "Local sample: Paired-End." >> "{log}"
            else
                touch {output.r2}
                echo "single" > {output.layout}
                echo "Local sample: Single-End." >> "{log}"
            fi
        fi
        """

# ------------------------------------------------------------------------------
# 4. Raw QC (Initial Check)
# ------------------------------------------------------------------------------

rule fastqc_raw:
    """Run FastQC on raw data (immediately after download/staging)."""
    container: config["docker_images"]["fastqc"]
    input:
        r1 = os.path.join(config["paths"]["data"], "{sample}/raw/reads_1.fastq.gz"),
        r2 = os.path.join(config["paths"]["data"], "{sample}/raw/reads_2.fastq.gz"),
        layout = os.path.join(config["paths"]["data"], "{sample}/raw/layout.txt")
    output:
        directory(os.path.join(config["paths"]["data"], "{sample}/fastqc_raw"))
    threads: config["resource_profiles"]["standard"]["threads"]
    resources:
        mem_mb = config["resource_profiles"]["standard"]["mem_mb"]
    log:
        os.path.join(config["paths"]["data"], "{sample}/logs/fastqc_raw.log")
    shell:
        """
        mkdir -p {output}
        layout=$(cat {input.layout})
        
        if [ "$layout" == "paired" ]; then
            fastqc --nogroup --threads {threads} -o {output} {input.r1} {input.r2} > "{log}" 2>&1
        else
            # Only run on R1, ignore the dummy R2
            fastqc --nogroup --threads {threads} -o {output} {input.r1} > "{log}" 2>&1
        fi
        """

rule multiqc_raw:
    """Aggregates FastQC reports from the Raw Data step."""
    container: config["docker_images"]["multiqc"]
    input:
        expand(os.path.join(config["paths"]["data"], "{sample}/fastqc_raw"), sample=SAMPLES)
    output:
        os.path.join(config["paths"]["data"], "MultiQC/multiqc_report.html")
    threads: config["resource_profiles"]["standard"]["threads"]
    resources:
        mem_mb = config["resource_profiles"]["standard"]["mem_mb"]
    log:
        "MultiQC.log"
    shell:
        """multiqc {input} -o $(dirname {output}) > "{log}" 2>&1"""