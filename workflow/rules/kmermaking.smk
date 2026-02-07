import os
import glob
import sys

# ==============================================================================
# Helper Functions for Checkpoints
# ==============================================================================

def get_all_chrom_names(wildcards):
    """
    Retrieves chromosome names after the genome split checkpoint.
    """
    # Wait for the checkpoint to finish
    checkpoint_output = checkpoints.step_04_split_masked_genome.get(**wildcards).output[0]
    
    # Find all fasta files in that directory
    search_path = os.path.join(checkpoint_output, "*.fasta")
    chrom_fastas = glob.glob(search_path)
    
    if not chrom_fastas:
        print(f"WARNING: No FASTA files found in {search_path}. Check permissions or binding.", file=sys.stderr)
        return []
        
    # Return basename without extension (e.g., "chr1" from "chr1.fasta")
    return [os.path.basename(f).replace(".fasta", "") for f in chrom_fastas]

def get_all_repeat_names(wildcards):
    """
    Retrieves repeat family names after the BED split checkpoint.
    """
    checkpoint_output = checkpoints.step_01_split_bed.get(**wildcards).output[0]
    
    search_path = os.path.join(checkpoint_output, "*.bed")
    bed_files = glob.glob(search_path)
    
    if not bed_files:
        print(f"WARNING: No BED files found in {search_path}.", file=sys.stderr)
        return []
        
    return [os.path.basename(f).replace(".bed", "") for f in bed_files]

# ==============================================================================
# Rules
# ==============================================================================

rule step_00_prepare_dirs:
    """Creates the base output directory for indexing."""
    output:
        directory(config["indexing"]["base_dir"])
    benchmark:
        "benchmarks/indexing/step_00_prepare_dirs.tsv"
    shell:
        "mkdir -p {output}"

# ------------------------------------------------------------------------------
# Phase 1: Preparation (BED Splitting & Fasta Extraction)
# ------------------------------------------------------------------------------

checkpoint step_01_split_bed:
    """
    Splits the master BED file and generates abundance statistics (ChrX/Y enrichment).
    """
    input:
        bed = config["resources"]["repeat_bed"]
    output:
        # The directory where split BED files go
        out_dir = directory(config["indexing"]["outputs"]["step_01_bed_split"]),
        # The directory where statistical summaries go (extracted from config or hardcoded)
        stats_dir = directory(config["indexing"]["outputs"]["step_01_bed_split"] + "_stats")
    benchmark:
        "benchmarks/indexing/step_01_split_bed.tsv"
    params:
        r_script = "workflow/scripts/split_bed_by_family.R"
    container:
        config["docker_images"]["r_scripting"]
    shell:
        """
        # Ensure output directories exist
        mkdir -p "{output.out_dir}"
        mkdir -p "{output.stats_dir}"

        # Call the R script
        # Note: We pass arguments matching the R script's 'argparse' definitions
        Rscript {params.r_script} \
            --input_bed "{input.bed}" \
            --output_dir "{output.out_dir}" \
            --intermediate "{output.stats_dir}"
        """

rule step_02_get_fasta_per_type:
    """
    Extracts FASTA sequences for each repeat family using bedtools.
    """
    input:
        genome = config["resources"]["genome_fasta"],
        bed = os.path.join(config["indexing"]["outputs"]["step_01_bed_split"], "{repeat}.bed")
    output:
        fasta = os.path.join(config["indexing"]["outputs"]["step_02_fasta_type"], "{repeat}.fasta")
    benchmark:
        "benchmarks/indexing/step_02_get_fasta_per_type_{repeat}.tsv"
    container:
        config["docker_images"]["bedtools"]
    shell:
        """
        bedtools getfasta -fi "{input.genome}" -bed "{input.bed}" -fo "{output.fasta}"
        """

# ------------------------------------------------------------------------------
# Phase 2: Genome Masking & Chromosome Splitting
# ------------------------------------------------------------------------------

rule step_03_mask_genome:
    """
    Creates a hard-masked version of the genome (Repeats masked with N).
    """
    input:
        genome = config["resources"]["genome_fasta"],
        bed = config["resources"]["repeat_bed"]
    output:
        masked_genome = config["indexing"]["outputs"]["step_03_masked_genome"]
    benchmark:
        "benchmarks/indexing/step_03_mask_genome.tsv"
    container:
        config["docker_images"]["bedtools"]
    shell:
        """
        bedtools maskfasta -fi "{input.genome}" -bed "{input.bed}" -fo "{output.masked_genome}"
        """

checkpoint step_04_split_masked_genome:
    """
    Splits the masked genome into individual chromosome FASTA files.
    """
    input:
        masked_genome = config["indexing"]["outputs"]["step_03_masked_genome"]
    output:
        directory(config["indexing"]["outputs"]["step_04_split_chrom"])
    benchmark:
        "benchmarks/indexing/step_04_split_masked_genome.tsv"
    container:
        config["docker_images"]["base_unix"] # Simple grep/awk
    shell:
        """
        mkdir -p {output}
        # Split fasta by header (simple implementation)
        awk '/^>/ {{OUT="{output}/" substr($0,2) ".fasta"}}; {{print >> OUT}}' {input.masked_genome}
        """

# ------------------------------------------------------------------------------
# Phase 3: K-mer Counting (Genome & Repeats)
# ------------------------------------------------------------------------------

rule step_05_count_chrom_kmers:
    """
    Counts k-mers in the non-repeat (masked) regions of each chromosome.
    """
    input:
        fasta = os.path.join(config["indexing"]["outputs"]["step_04_split_chrom"], "{chrom}.fasta")
    output:
        txt = os.path.join(config["indexing"]["outputs"]["step_05_chrom_counts"], "{chrom}.txt")
    benchmark:
        "benchmarks/indexing/step_05_jellyfish_{chrom}.tsv"
    params:
        k = config["indexing"]["kmer_len"],
        size = config["indexing"]["hash_size"]
    threads: 
        config["resource_profiles"]["heavy"]["threads"]
    resources:
        mem_mb = config["resource_profiles"]["heavy"]["mem_mb"],
        tmpdir = os.path.join(config["paths"]["data"], "tmp")
    container:
        config["docker_images"]["jellyfish"]
    shell:
        """
        # 1. Create a temp file safely inside the job directory
        TMP_DB=$(mktemp -u --tmpdir={resources.tmpdir} jf_XXXXXX)
        # 2. Count with --canonical
        jellyfish count \
            --mer-len {params.k} \
            --size {params.size} \
            --threads {threads} \
            --canonical \
            --output $TMP_DB \
            {input.fasta}
        jellyfish dump --column $TMP_DB > "{output.txt}"
        rm -f $TMP_DB
        """

rule step_06_dump_genome_kmers:
    """
    Dumps all non-repeat k-mers into a single text file.
    This aggregates results from ALL chromosomes.
    """
    input:
        # Aggregation using the checkpoint function
        lambda w: expand(
            os.path.join(config["indexing"]["outputs"]["step_05_chrom_counts"], "{chrom}.txt"),
            chrom=get_all_chrom_names(w)
        )
    output:
        txt = config["indexing"]["outputs"]["step_06_genome_kmers"]
    benchmark:
        "benchmarks/indexing/step_06_dump_genome_kmers.tsv"
    params:
        input_dir = config["indexing"]["outputs"]["step_05_chrom_counts"]
    threads: 
        config["resource_profiles"]["heavy"]["threads"]
    resources:
        mem_mb = config["resource_profiles"]["heavy"]["mem_mb"]
    container:
        config["docker_images"]["base_unix"]
    shell:
        """
        find "{params.input_dir}" -name "*.txt" -print0 | \\
        xargs -0 cat | \\
        awk '{{print $1}}' | \\
        LC_ALL=C sort -u --parallel={threads} -S {resources.mem_mb}M > {output.txt}
        """

rule step_07_create_genome_db:
    """
    Loads non-repeat k-mers into DuckDB for fast querying.
    """
    input:
        txt = config["indexing"]["outputs"]["step_06_genome_kmers"]
    output:
        db = config["indexing"]["outputs"]["step_07_genome_db"]
    benchmark:
        "benchmarks/indexing/step_07_duckdb.tsv"
    params:
        sanitized_path = lambda w, input: input.txt.replace("'", "''")
    resources:
        mem_mb = config["resource_profiles"]["heavy"]["mem_mb"]
    container:
        config["docker_images"]["duckdb"]
    shell:
        """
        duckdb :memory: "COPY ( \
            SELECT column0 AS kmer \
            FROM read_csv('{params.sanitized_path}', header=False, columns={{'column0':'VARCHAR'}}) \
        ) TO '{output}' (FORMAT PARQUET, COMPRESSION 'SNAPPY')"
        """

rule step_08_count_repeat_kmers:
    """
    Counts k-mers for each repeat family type.
    """
    input:
        fasta = os.path.join(config["indexing"]["outputs"]["step_02_fasta_type"], "{repeat}.fasta")
    output:
        txt = os.path.join(config["indexing"]["outputs"]["step_08_repeat_counts"], "{repeat}.txt")
    benchmark:
        "benchmarks/indexing/step_08_count_repeat_kmers_{repeat}.tsv"
    params:
        k = config["indexing"]["kmer_len"],
        size = config["indexing"]["hash_size"]
    threads:
        config["resource_profiles"]["standard"]["threads"]
    resources:
        mem_mb = config["resource_profiles"]["standard"]["mem_mb"],
        tmpdir = os.path.join(config["paths"]["data"], "tmp")
    container:
        config["docker_images"]["jellyfish"]
    shell:
        """
        TMP_DB=$(mktemp -u --tmpdir={resources.tmpdir} jf_rep_XXXXXX)
        
        jellyfish count \
            --mer-len {params.k} \
            --size {params.size} \
            --threads {threads} \
            --canonical \
            --lower-count=1 \
            --output $TMP_DB \
            "{input.fasta}"

        jellyfish dump --column $TMP_DB > "{output.txt}"
        
        rm -f $TMP_DB
        """

rule step_09_create_global_repeat_db:
    """
    Aggregates ALL repeat k-mers into a single DuckDB/Parquet file.
    This enables checking "Is this k-mer present in ANY other repeat?"
    """
    input:
        lambda w: expand(
            os.path.join(config["indexing"]["outputs"]["step_08_repeat_counts"], "{repeat}.txt"),
            repeat=get_all_repeat_names(w)
        )
    output:
        db = config["indexing"]["outputs"]["step_09_global_db"]
    benchmark:
        "benchmarks/indexing/step_09_create_global_repeat_db.tsv"
    params:
        input_dir = config["indexing"]["outputs"]["step_08_repeat_counts"],
        temp_txt = "temp_global_sum.txt"
    threads:
        config["resource_profiles"]["heavy"]["threads"]
    resources:
        mem_mb = config["resource_profiles"]["heavy"]["mem_mb"]
    container:
        config["docker_images"]["duckdb"]
    shell:
        """
        export LC_ALL=C
        
        find {params.input_dir} -name "*.txt" -print0 | \
        xargs -0 cat | \
        sort --parallel={threads} -S {resources.mem_mb}M -k1,1 | \
        awk 'BEGIN {{ prev=""; sum=0; }} {{ if ($1 == prev) {{ sum += $2; }} else {{ if (prev != "") print prev, sum; prev = $1; sum = $2; }} }} END {{ if (prev != "") print prev, sum; }}' \
        > {params.temp_txt}
        
        # 2. Load into DuckDB
        # Schema: kmer (VARCHAR), total_count (INTEGER)
        duckdb :memory: "COPY ( \
            SELECT column0 AS kmer, column1 AS total_count \
            FROM read_csv('{params.temp_txt}', delim=' ', header=False, columns={{'column0':'VARCHAR', 'column1':'INTEGER'}}) \
        ) TO '{output.db}' (FORMAT PARQUET, COMPRESSION 'SNAPPY')"
        
        rm -f {params.temp_txt}
        """

rule step_10_identify_unique_repeat_kmers:
    """
    Filters k-mers that are unique (or nearly unique) to a specific repeat family.
    """
    input:
        txt = os.path.join(config["indexing"]["outputs"]["step_08_repeat_counts"], "{repeat}.txt"),
        global_db = config["indexing"]["outputs"]["step_09_global_db"]
    output:
        unique_txt = os.path.join(config["indexing"]["outputs"]["step_10_unique_repeat"], "{repeat}.txt")
    benchmark:
        "benchmarks/indexing/step_10_identify_unique_repeat_kmers_{repeat}.tsv"
    params:
        sanitized_path = lambda w, input: input.txt.replace("'", "''")
    threads:
        config["resource_profiles"]["light"]["threads"]
    resources:
        mem_mb = config["resource_profiles"]["light"]["mem_mb"]
    container:
        config["docker_images"]["duckdb"]
    shell:
        """      
        duckdb :memory: "SET memory_limit='{resources.mem_mb}MB'; \
            SET threads={threads}; \
            COPY ( \
            WITH local_counts AS ( \
                SELECT column0 as kmer, column1 as count \
                FROM read_csv('{params.sanitized_path}', delim=' ', header=False, columns={{'column0':'VARCHAR', 'column1':'INTEGER'}}) \
            ) \
            SELECT lc.kmer \
            FROM local_counts lc \
            INNER JOIN read_parquet('{input.global_db}') gc ON lc.kmer = gc.kmer \
            WHERE ((gc.total_count - lc.count)::DOUBLE / gc.total_count) < 0.01 \
        ) TO '{output.unique_txt}' (FORMAT CSV, DELIMITER '\\t', HEADER TRUE)"
        """

rule step_11_remove_genomic_kmers:
    """
    Removes k-mers that appear in the non-repeat genome (False Positives).
    """
    input:
        repeat_txt = os.path.join(config["indexing"]["outputs"]["step_10_unique_repeat"], "{repeat}.txt"),
        genome_db = config["indexing"]["outputs"]["step_07_genome_db"]
    output:
        clean_fasta = os.path.join(config["indexing"]["outputs"]["step_11_genomic_removed"], "{repeat}.fasta")
    benchmark:
        "benchmarks/indexing/step_11_remove_genomic_kmers_{repeat}.tsv"
    params:
        sanitized_path = lambda w, input: input.repeat_txt.replace("'", "''")
    threads:
        config["resource_profiles"]["light"]["threads"]
    resources:
        mem_mb = config["resource_profiles"]["light"]["mem_mb"]
    container:
        config["docker_images"]["duckdb"]
    shell:
        """
        duckdb :memory: "SET memory_limit='{resources.mem_mb}MB'; \
            SET threads={threads}; \
            COPY ( \
            SELECT local.kmer \
            FROM read_csv('{params.sanitized_path}', delim='\\t', header=TRUE) AS local \
            LEFT JOIN read_parquet('{input.genome_db}') AS global ON local.kmer = global.kmer \
            WHERE global.kmer IS NULL \
        ) TO '/dev/stdout' (FORMAT CSV, HEADER FALSE)" \
        | awk -v name="{wildcards.repeat}" '{{print ">" name "_" NR "\\n" $0}}' > "{output.clean_fasta}"
        """

# ------------------------------------------------------------------------------
# Phase 4: Transcriptome Filtering
# ------------------------------------------------------------------------------

rule step_12_generate_transcriptome_fasta:
    """
    Generates a standardized transcriptome FASTA from Genome + GTF.
    This ensures headers match the GTF (e.g., ENST IDs) for downstream mapping.
    """
    input:
        genome = config["resources"]["genome_fasta"],
        gtf = config["resources"]["genome_gtf"]
    output:
        fasta = config["indexing"]["generated_transcriptome"]
    benchmark:
        "benchmarks/indexing/step_12_generate_transcriptome_fasta.tsv"
    container: config["docker_images"]["gffread"]
    shell:
        """
        # -w outputs the transcript sequences
        # -g provides the genomic reference for coordinate lookup
        gffread -w "{output.fasta}" -g "{input.genome}" "{input.gtf}"
        """

rule step_12_process_transcripts:
    """
    Counts k-mers in the canonical transcriptome to remove further FPs.
    """
    input:
        transcripts = config["indexing"]["generated_transcriptome"]
    output:
        txt = config["indexing"]["outputs"]["step_12_trans_txt"]
    benchmark:
        "benchmarks/indexing/step_12_process_transcripts.tsv"
    params:
        k = config["indexing"]["kmer_len"],
        size = config["indexing"]["hash_size"]
    threads:
        config["resource_profiles"]["heavy"]["threads"]
    resources:
        mem_mb = config["resource_profiles"]["heavy"]["mem_mb"],
        tmpdir = os.path.join(config["paths"]["data"], "tmp")
    container:
        config["docker_images"]["jellyfish"]
    shell:
        """
        TMP_DB=$(mktemp -u --tmpdir={resources.tmpdir} jf_trans_XXXXXX)
        
        jellyfish count \
            --mer-len {params.k} \
            --size {params.size} \
            --threads {threads} \
            --canonical \
            --output $TMP_DB \
            "{input.transcripts}"
            
        jellyfish dump --column $TMP_DB > "{output.txt}"
        
        rm -f $TMP_DB
        """

rule step_13_create_transcript_db:
    input:
        txt = config["indexing"]["outputs"]["step_12_trans_txt"]
    output:
        db = config["indexing"]["outputs"]["step_13_trans_db"]
    benchmark:
        "benchmarks/indexing/step_13_create_transcript_db.tsv"
    params:
        sanitized_path = lambda w, input: input.txt.replace("'", "''")
    resources:
        mem_mb = config["resource_profiles"]["heavy"]["mem_mb"]
    container:
        config["docker_images"]["duckdb"]
    shell:
        """
        duckdb :memory: "COPY ( \
            SELECT column0 AS kmer \
            FROM read_csv('{params.sanitized_path}', header=False, delim=' ', columns={{'column0':'VARCHAR', 'column1':'INTEGER'}}) \
        ) TO '{output.db}' (FORMAT PARQUET, COMPRESSION 'SNAPPY')"
        """

rule step_14_final_filtering:
    """
    Final Filter: Remove k-mers present in transcripts.
    ADAPTED ORIGINAL: Preserves original FASTA headers using 'paste - -'.
    """
    input:
        clean_fasta = os.path.join(config["indexing"]["outputs"]["step_11_genomic_removed"], "{repeat}.fasta"),
        trans_db = config["indexing"]["outputs"]["step_13_trans_db"]
    output:
        final_fasta = os.path.join(config["indexing"]["outputs"]["step_14_final_fasta"], "{repeat}.fasta")
    benchmark:
        "benchmarks/indexing/step_14_final_filtering_{repeat}.tsv"
    threads:
        config["resource_profiles"]["light"]["threads"]
    resources:
        mem_mb = config["resource_profiles"]["light"]["mem_mb"]
    container:
        config["docker_images"]["duckdb"]
    shell:
        """
        if [ -s "{input.clean_fasta}" ]; then
            paste - - < "{input.clean_fasta}" \
            | duckdb :memory: "SET memory_limit='{resources.mem_mb}MB'; \
                SET threads={threads}; \
                COPY ( \
                SELECT local.column0, local.column1 \
                FROM read_csv('/dev/stdin', delim='\\t', header=FALSE, columns={{'column0':'VARCHAR', 'column1':'VARCHAR'}}) AS local \
                LEFT JOIN read_parquet('{input.trans_db}') AS global ON local.column1 = global.kmer \
                WHERE global.kmer IS NULL \
            ) TO '/dev/stdout' (FORMAT CSV, DELIMITER '\\t', HEADER FALSE)" \
            | tr '\\t' '\\n' > "{output.final_fasta}"
        else
            # Handle empty input safely
            touch "{output.final_fasta}"
        fi
        """

# ------------------------------------------------------------------------------
# Phase 6: Reporting & Cleanup
# ------------------------------------------------------------------------------

rule step_15_generate_summary:
    """
    Generates a summary table showing how many k-mers survived filtering.
    """
    input:
        files_before = lambda w: expand(
            os.path.join(config["indexing"]["outputs"]["step_11_genomic_removed"], "{repeat}.fasta"),
            repeat=get_all_repeat_names(w)
        ),
        files_after = lambda w: expand(
            os.path.join(config["indexing"]["outputs"]["step_14_final_fasta"], "{repeat}.fasta"),
            repeat=get_all_repeat_names(w)
        )
    output:
        summary = config["indexing"]["outputs"]["step_15_summary"]
    benchmark:
        "benchmarks/indexing/step_15_generate_summary.tsv"
    params:
        dir_before = config["indexing"]["outputs"]["step_11_genomic_removed"],
        dir_after = config["indexing"]["outputs"]["step_14_final_fasta"]
    container:
        config["docker_images"]["base_unix"]
    shell:
        """
        # Header: Name, After (Final), Before (Transcript Filter)
        printf "Name\\tAfter\\tBefore\\n" > "{output.summary}"
        
        for file_path_A in "{params.dir_before}"/*; do
            filename=$(basename "$file_path_A")
            file_path_B="{params.dir_after}/$filename"
            
            if [ -f "$file_path_B" ]; then
                count_before=$(grep -c "^>" "$file_path_A" || true)
                count_after=$(grep -c "^>" "$file_path_B" || true)
                name="${{filename%.*}}"
                printf "%s\\t%s\\t%s\\n" "$name" "$count_after" "$count_before" >> "{output.summary}"
            fi
        done
        """

rule step_16_cleanup:
    """
    Removes empty FASTA files from the final output directory.
    Run only after summary is generated.
    """
    input:
        config["indexing"]["outputs"]["step_15_summary"]
    output:
        touch(config["indexing"]["outputs"]["step_16_cleanup_done"])
    benchmark:
        "benchmarks/indexing/step_16_cleanup.tsv"
    params:
        final_dir = config["indexing"]["outputs"]["step_14_final_fasta"]
    shell:
        """
        find {params.final_dir} -name "*.fasta" -size 0 -delete
        """
