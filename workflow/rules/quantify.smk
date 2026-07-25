import os
import json
import pandas as pd

# ------------------------------------------------------------------------------
# 1. Trimming & QC (Fastp)
# ------------------------------------------------------------------------------

rule fastp:
    """
    Trim adapters and low quality bases using fastp.
    Input comes from download.smk (normalized raw/reads_*.fastq.gz).
    """
    container: config["docker_images"]["fastp"]
    input:
        r1 = os.path.join(config["paths"]["data"], "{sample}/raw/reads_1.fastq.gz"),
        r2 = os.path.join(config["paths"]["data"], "{sample}/raw/reads_2.fastq.gz"),
        layout = os.path.join(config["paths"]["data"], "{sample}/raw/layout.txt")
    output:
        r1 = os.path.join(config["paths"]["data"], "{sample}/fastp/reads_1.trimmed.fastq.gz"),
        r2 = os.path.join(config["paths"]["data"], "{sample}/fastp/reads_2.trimmed.fastq.gz"),
        json = os.path.join(config["paths"]["data"], "{sample}/fastp/report.json"),
        html = os.path.join(config["paths"]["data"], "{sample}/fastp/report.html"),
        done = touch(os.path.join(config["paths"]["data"], "{sample}/fastp/trimming.done"))
    log:
        os.path.join(config["paths"]["data"], "{sample}/logs/fastp.log")
    threads: config["resource_profiles"]["standard"]["threads"]
    resources:
        mem_mb = config["resource_profiles"]["standard"]["mem_mb"]
    shell:
        """
        layout=$(cat {input.layout})
        
        # Ensure output directory exists
        mkdir -p $(dirname {output.r1})

        if [ "$layout" == "paired" ]; then
            echo "Running Fastp in Paired-End mode..." > "{log}"
            fastp \
                -i {input.r1} -I {input.r2} \
                -o {output.r1} -O {output.r2} \
                -j {output.json} -h {output.html} \
                --thread {threads} \
                --detect_adapter_for_pe \
                2>> "{log}"
        else
            echo "Running Fastp in Single-End mode..." > "{log}"
            fastp \
                -i {input.r1} \
                -o {output.r1} \
                -j {output.json} -h {output.html} \
                --thread {threads} \
                2>> "{log}"
                
            # Create a dummy R2 file so Snakemake job satisfies output requirements
            # We touch it to verify it exists but is empty
            touch {output.r2}
        fi
        """

rule fastqc_trimmed:
    """Run FastQC on trimmed data to verify quality improvements."""
    container: config["docker_images"]["fastqc"]
    input:
        r1 = os.path.join(config["paths"]["data"], "{sample}/fastp/reads_1.trimmed.fastq.gz"),
        r2 = os.path.join(config["paths"]["data"], "{sample}/fastp/reads_2.trimmed.fastq.gz"),
        layout_flag = os.path.join(config["paths"]["data"], "{sample}/fastp/trimming.done"),
    output:
        directory(os.path.join(config["paths"]["data"], "{sample}/fastp-fastqc"))
    threads: config["resource_profiles"]["standard"]["threads"]
    resources:
        mem_mb = config["resource_profiles"]["standard"]["mem_mb"]
    log:
        os.path.join(config["paths"]["data"], "{sample}/logs/fastqc_trimmed.log")
    shell:
        """
        mkdir -p {output}
        
        # Always QC R1
        fastqc --nogroup --threads {threads} -o {output} {input.r1} > "{log}" 2>&1
        
        # Check if R2 is non-empty (it is 0 bytes if SE)
        if [ -s {input.r2} ]; then
             fastqc --nogroup --threads {threads} -o {output} {input.r2} >> "{log}" 2>&1
        fi
        """

rule multiqc_trimmed:
    """Aggregates FastQC reports from the Trimmed Data step."""
    container: config["docker_images"]["multiqc"]
    input:
        expand(os.path.join(config["paths"]["data"], "{sample}/fastp-fastqc"), sample=SAMPLES)
    output:
        os.path.join(config["paths"]["data"], "fastp-MultiQC/multiqc_report.html")
    threads: config["resource_profiles"]["standard"]["threads"]
    resources:
        mem_mb = config["resource_profiles"]["standard"]["mem_mb"]
    log:
        "logs/fastp-MultiQC.log"
    shell:
        """multiqc {input} -o $(dirname {output}) > "{log}" 2>&1"""

# ------------------------------------------------------------------------------
# 2. Repeat Quantification (Rust Tool)
# ------------------------------------------------------------------------------

rule quantify_repeats:
    """
    Quantifies repetitive elements with optional 1:1 read assignment tracking.
    """
    container: config["docker_images"]["repeatquant"]
    input:
        r1 = os.path.join(config["paths"]["data"], "{sample}/fastp/reads_1.trimmed.fastq.gz"),
        r2 = os.path.join(config["paths"]["data"], "{sample}/fastp/reads_2.trimmed.fastq.gz"),
        layout = os.path.join(config["paths"]["data"], "{sample}/raw/layout.txt"),
        index = config["quantification"]["repeat_index_path"]
    output:
        counts = os.path.join(config["paths"]["data"], "{sample}/quant_repeats/counts.tsv"),
        assignments = os.path.join(config["paths"]["data"], "{sample}/quant_repeats/assignments.txt") if config["quantification"].get("export_assignments", False) else []
    params:
        min_matches = config["quantification"].get("min_kmer_matches", 1),
        assign_flag = lambda wildcards, output: f"--assignments-out {output.assignments}" if hasattr(output, "assignments") and output.assignments else ""
    threads: config["resource_profiles"]["extreme"]["threads"]
    benchmark: os.path.join(config["paths"]["benchmarks"], "quant_repeats","{sample}.txt")
    resources:
        mem_mb = config["resource_profiles"]["extreme"]["mem_mb"]
    shell:
        """
        layout=$(cat {input.layout})

        if [ "$layout" == "paired" ]; then
            repeat_counter quant \
                --index {input.index} \
                --r1 {input.r1} \
                --r2 {input.r2} \
                --threads {threads} \
                --min-matches {params.min_matches} \
                {params.assign_flag} \
                --output {output.counts}
        else
            repeat_counter quant \
                --index {input.index} \
                --r1 {input.r1} \
                --threads {threads} \
                --min-matches {params.min_matches} \
                {params.assign_flag} \
                --output {output.counts}
        fi
        """

rule kmer_aggregate_matrix:
    """
    Aggregates Rust RepeatQuant outputs into Raw and CPM matrices.
    """
    #container: config["docker_images"]["pand"]
    input:
        cols = expand(os.path.join(config["paths"]["data"], "{sample}/quant_repeats/counts.tsv"), sample=SAMPLES),
        json_reports = expand(os.path.join(config["paths"]["data"], "{sample}/fastp/report.json"), sample=SAMPLES)
    output:
        matrix_cpm = os.path.join(config["paths"]["results"], "kmer_counts_cpm.tsv"),
        matrix_raw = os.path.join(config["paths"]["results"], "kmer_counts_raw.tsv"),
        stats = os.path.join(config["paths"]["results"], "sample_total_reads.tsv")
    log:
        "logs/kmer_aggregate.log"
    run:
        dfs_cpm = []
        dfs_raw = []
        read_stats = []
        json_path_map = {}
        
        # Map sample names to fastp logs
        # Logic: input is data/{sample}/fastp/report.json
        for p in input.json_reports:
            # dirname(p) -> data/{sample}/fastp
            # dirname(dirname(p)) -> data/{sample}
            # basename(...) -> {sample}
            name_key = os.path.basename(os.path.dirname(os.path.dirname(p)))
            json_path_map[name_key] = p
        
        for col_file in input.cols:
            # Logic: input is data/{sample}/quant_repeats/counts.tsv
            name_wildcard = os.path.basename(os.path.dirname(os.path.dirname(col_file)))
            
            # Defensive check
            if name_wildcard not in json_path_map:
                print(f"Warning: No JSON found for {name_wildcard}")
                continue

            json_path = json_path_map[name_wildcard]
            
            # 1. Get Library Size (Total Reads)
            with open(json_path, 'r') as f_json:
                fastp_data = json.load(f_json)
            
            # Use 'total_reads' from 'after_filtering'
            total_reads = float(fastp_data['summary']['after_filtering']['total_reads'])

            read_stats.append({"sample_name": name_wildcard, "total_reads": int(total_reads)})
    
            # 2. Read Rust Output
            # Format: RepeatType \t Count
            df = pd.read_csv(col_file, sep="\t")
            if "RepeatType" not in df.columns:
                 # Handle cases where header might be missing or different
                 # Assuming 1st col is index, 2nd is count if headers fail
                 pass 
            
            df = df.set_index("RepeatType")
            
            # 3. Create Raw Count Series
            df_raw = df.rename(columns={"Count": name_wildcard})
            dfs_raw.append(df_raw)

            # 4. Create CPM Series
            # CPM = (Count / Total_Reads) * 1,000,000
            df_cpm = df_raw.copy()
            df_cpm[name_wildcard] = (df_cpm[name_wildcard] / total_reads) * 1_000_000
            dfs_cpm.append(df_cpm)

        # 5. Merge and Save RAW
        if dfs_raw:
            merged_raw = pd.concat(dfs_raw, axis=1).fillna(0.0)
            merged_raw.to_csv(output.matrix_raw, sep="\t", index_label="RepeatType")
            
            # 6. Merge and Save CPM
            merged_cpm = pd.concat(dfs_cpm, axis=1).fillna(0.0)
            merged_cpm.to_csv(output.matrix_cpm, sep="\t", index_label="RepeatType")

            df_stats = pd.DataFrame(read_stats)
            df_stats = df_stats.sort_values("sample_name") 
            df_stats.to_csv(output.stats, sep="\t", index=False)

# ------------------------------------------------------------------------------
# 3. Gene Quantification (Kallisto)
# ------------------------------------------------------------------------------

rule generate_tx2gene:
    """
    Creates a Transcript-to-Gene mapping file from the GTF.
    This allows aggregating Kallisto transcript counts into gene counts (Gene Symbol).
    """
    container: config["docker_images"]["gawk"]
    input:
        gtf = config["resources"]["genome_gtf"]
    output:
        map = config["quantification"]["tx2gene"]
    log:
        "logs/generate_tx2gene.log"
    shell:
        r"""
        awk 'BEGIN {{FS="\t"; OFS="\t"}} \
            $3 == "transcript" {{ \
                tx_id = ""; \
                gn_id = ""; \
                n = split($9, fields, ";"); \
                for (i=1; i<=n; i++) {{ \
                    gsub(/^[ \t\xA0]+|[ \t\xA0]+$/, "", fields[i]); \
                    if (fields[i] ~ /^transcript_id/) {{ \
                        match(fields[i], /"([^"]+)"/, tx_val); \
                        tx_id = tx_val[1]; \
                    }} \
                    if (fields[i] ~ /^gene_name/) {{ \
                        match(fields[i], /"([^"]+)"/, gn_val); \
                        gn_id = gn_val[1]; \
                    }} \
                }} \
                if (tx_id != "" && gn_id != "") {{ \
                    print tx_id, gn_id \
                }} \
            }}' {input.gtf} 2> "{log}" | sort -u > {output.map}
        """

rule kallisto_index:
    """
    Index the transcriptome for Kallisto.
    Uses the EXISTING transcriptome FASTA as kmer filtering step.
    """
    container: config["docker_images"]["kallisto"]
    input:
        transcriptome = config["indexing"]["generated_transcriptome"]
    output:
        index = "resources/kallisto.idx"
    threads: config["resource_profiles"]["heavy"]["threads"]
    resources:
        mem_mb = config["resource_profiles"]["heavy"]["mem_mb"]
    log:
        "logs/kallisto_index.log"
    shell:
        """
        kallisto index -i {output.index} {input.transcriptome} > "{log}" 2>&1
        """

rule kallisto_quant:
    """Quantify gene expression using Kallisto."""
    container: config["docker_images"]["kallisto"]
    input:
        r1 = os.path.join(config["paths"]["data"], "{sample}/fastp/reads_1.trimmed.fastq.gz"),
        r2 = os.path.join(config["paths"]["data"], "{sample}/fastp/reads_2.trimmed.fastq.gz"),
        layout = os.path.join(config["paths"]["data"], "{sample}/raw/layout.txt"),
        index = "resources/kallisto.idx"
    output:
        # Kallisto outputs 'abundance.tsv' and 'abundance.h5' in a directory
        outdir = directory(os.path.join(config["paths"]["data"], "{sample}/quant_genes")),
        tsv = os.path.join(config["paths"]["data"], "{sample}/quant_genes/abundance.tsv")
    params:
        # Defaults for Single End if not estimated
        frag_len = config["resources"]["kallisto"].get("fragment_len", 200),
        frag_sd = config["resources"]["kallisto"].get("fragment_sd", 20)
    threads: config["resource_profiles"]["heavy"]["threads"]
    resources:
        mem_mb = config["resource_profiles"]["heavy"]["mem_mb"]
    log:
        os.path.join(config["paths"]["data"], "{sample}/logs/kallisto.log")
    shell:
        """
        layout=$(cat {input.layout})
        
        if [ "$layout" == "paired" ]; then
            kallisto quant \
                -i {input.index} \
                -o {output.outdir} \
                --threads {threads} \
                {input.r1} {input.r2} > "{log}" 2>&1
        else
            # SE requires --single -l -s
            kallisto quant \
                -i {input.index} \
                -o {output.outdir} \
                --threads {threads} \
                --single -l {params.frag_len} -s {params.frag_sd} \
                {input.r1} > "{log}" 2>&1
        fi

        rm -f {output.outdir}/abundance.h5
        """

rule aggr_kallisto_genes:
    """
    Aggregates Kallisto outputs using tximport (R).
    Uses the dynamically generated tx2gene mapping.
    """
    container: config["docker_images"]["r_tximport"]
    input:
        quant_files = expand(os.path.join(config["paths"]["data"], "{sample}/quant_genes/abundance.tsv"), sample=SAMPLES),
        tx2gene = config["quantification"]["tx2gene"]
    output:
        counts = os.path.join(config["paths"]["results"], "gene_counts_raw.tsv"),
        tpm = os.path.join(config["paths"]["results"], "gene_counts_tpm.tsv")
    params:
        sample_names = SAMPLES
    log:
        "logs/tximport.log"
    script:
        "../scripts/tximport.R"
