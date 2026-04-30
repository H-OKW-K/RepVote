import os
import glob
import sys

# ==============================================================================
# Helper Functions for Checkpoints
# ==============================================================================

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
    shell:
        "mkdir -p {output}"

# ------------------------------------------------------------------------------
# Phase 1: Preparation
# ------------------------------------------------------------------------------

checkpoint step_01_split_bed:
    """Splits the master BED file."""
    input:
        bed = config["resources"]["repeat_bed"]
    output:
        out_dir = directory(config["indexing"]["outputs"]["step_01_bed_split"]),
        stats_dir = directory(config["indexing"]["outputs"]["step_01_bed_split"] + "_stats")
    params:
        r_script = "workflow/scripts/split_bed_by_family.R"
    container:
        config["docker_images"]["r_scripting"]
    shell:
        """
        mkdir -p "{output.out_dir}" "{output.stats_dir}"
        Rscript {params.r_script} \
            --input_bed "{input.bed}" \
            --output_dir "{output.out_dir}" \
            --intermediate "{output.stats_dir}"
        """

rule step_02_get_fasta_per_type:
    """Extracts FASTA sequences for each repeat family."""
    input:
        genome = config["resources"]["genome_fasta"],
        bed = os.path.join(config["indexing"]["outputs"]["step_01_bed_split"], "{repeat}.bed")
    output:
        fasta = os.path.join(config["indexing"]["outputs"]["step_02_fasta_type"], "{repeat}.fasta")
    container:
        config["docker_images"]["bedtools"]
    shell:
        "bedtools getfasta -fi '{input.genome}' -bed '{input.bed}' -fo '{output.fasta}'"

rule step_03_mask_genome:
    """Creates a hard-masked version of the genome."""
    input:
        genome = config["resources"]["genome_fasta"],
        bed = config["resources"]["repeat_bed"]
    output:
        masked_genome = config["indexing"]["outputs"]["step_03_masked_genome"]
    container:
        config["docker_images"]["bedtools"]
    shell:
        "bedtools maskfasta -fi '{input.genome}' -bed '{input.bed}' -fo '{output.masked_genome}'"

rule step_04_generate_transcriptome:
    """Generates standard transcriptome FASTA."""
    input:
        genome = config["resources"]["genome_fasta"],
        gtf = config["resources"]["genome_gtf"]
    output:
        fasta = config["indexing"]["generated_transcriptome"]
    container: 
        config["docker_images"]["gffread"]
    shell:
        "gffread -w '{output.fasta}' -g '{input.genome}' '{input.gtf}'"

rule step_05_get_all_repeats_fasta:
    """Extracts ALL repeats into a single FASTA for global counting."""
    input:
        genome = config["resources"]["genome_fasta"],
        bed = config["resources"]["repeat_bed"]
    output:
        fasta = config["indexing"]["outputs"]["step_05_all_repeats_fasta"]
    container:
        config["docker_images"]["bedtools"]
    shell:
        "bedtools getfasta -fi '{input.genome}' -bed '{input.bed}' -fo '{output.fasta}'"

# ------------------------------------------------------------------------------
# Phase 2: Core KMC3 Indexing (Genomes, Transcripts, Global Repeats)
# ------------------------------------------------------------------------------

rule step_06a_build_genome_kmc_db:
    """Builds foundational KMC database for the genome."""
    input:
        fasta = config["indexing"]["outputs"]["step_03_masked_genome"]
    output:
        pre = config["indexing"]["outputs"]["genome_db"] + ".kmc_pre",
        suf = config["indexing"]["outputs"]["genome_db"] + ".kmc_suf"
    params:
        k = config["indexing"]["kmer_len"],
        prefix = config["indexing"]["outputs"]["genome_db"]
    threads: config["resource_profiles"]["heavy"]["threads"]
    resources:
        tmpdir = os.path.join(config["paths"]["data"], "tmp_genome")
    container: config["docker_images"]["kmc"]
    shell:
        """
        mkdir -p {resources.tmpdir}
        kmc -k{params.k} -ci1 -fm -t{threads} "{input.fasta}" "{params.prefix}" {resources.tmpdir}
        """

rule step_06b_build_trans_kmc_db:
    """Builds foundational KMC database for the transcriptome."""
    input:
        fasta = config["indexing"]["generated_transcriptome"]
    output:
        pre = config["indexing"]["outputs"]["trans_db"] + ".kmc_pre",
        suf = config["indexing"]["outputs"]["trans_db"] + ".kmc_suf"
    params:
        k = config["indexing"]["kmer_len"],
        prefix = config["indexing"]["outputs"]["trans_db"]
    threads: config["resource_profiles"]["heavy"]["threads"]
    resources:
        tmpdir = os.path.join(config["paths"]["data"], "tmp_trans")
    container: config["docker_images"]["kmc"]
    shell:
        """
        mkdir -p {resources.tmpdir}
        kmc -k{params.k} -ci1 -fm -t{threads} "{input.fasta}" "{params.prefix}" {resources.tmpdir}
        """

rule step_06c_build_rep_kmc_db:
    """Builds foundational KMC database for all repeats."""
    input:
        fasta = config["indexing"]["outputs"]["step_05_all_repeats_fasta"]
    output:
        pre = config["indexing"]["outputs"]["all_rep_db"] + ".kmc_pre",
        suf = config["indexing"]["outputs"]["all_rep_db"] + ".kmc_suf"
    params:
        k = config["indexing"]["kmer_len"],
        prefix = config["indexing"]["outputs"]["all_rep_db"]
    threads: config["resource_profiles"]["heavy"]["threads"]
    resources:
        tmpdir = os.path.join(config["paths"]["data"], "tmp_rep")
    container: config["docker_images"]["kmc"]
    shell:
        """
        mkdir -p {resources.tmpdir}
        kmc -k{params.k} -ci1 -fm -t{threads} "{input.fasta}" "{params.prefix}" {resources.tmpdir}
        """

rule step_06d_build_global_valid_repeats:
    """Filters all repeats against the genome/transcriptome EXACTLY ONCE."""
    input:
        gen_pre = config["indexing"]["outputs"]["genome_db"] + ".kmc_pre",
        trans_pre = config["indexing"]["outputs"]["trans_db"] + ".kmc_pre",
        rep_pre = config["indexing"]["outputs"]["all_rep_db"] + ".kmc_pre"
    output:
        # A master text dictionary of every allowed k-mer and its GLOBAL count
        valid_txt = config["indexing"]["outputs"]["all_rep_db"] + "_valid.txt"
    params:
        gen_db = config["indexing"]["outputs"]["genome_db"],
        trans_db = config["indexing"]["outputs"]["trans_db"],
        rep_db = config["indexing"]["outputs"]["all_rep_db"],
        temp_db = lambda w, resources: f"{resources.tmpdir}/global_temp",
        valid_db = lambda w, resources: f"{resources.tmpdir}/global_valid"
    threads: config["resource_profiles"]["heavy"]["threads"]
    resources:
        tmpdir = os.path.join(config["paths"]["data"], "tmp_global")
    container: config["docker_images"]["kmc"]
    shell:
        """
        mkdir -p {resources.tmpdir}
        
        # 1. Pairwise subtract the Genome from All Repeats
        kmc_tools -t{threads} simple "{params.rep_db}" "{params.gen_db}" kmers_subtract "{params.temp_db}"
        
        # 2. Pairwise subtract the Transcriptome from the surviving repeats
        kmc_tools -t{threads} simple "{params.temp_db}" "{params.trans_db}" kmers_subtract "{params.valid_db}"
        
        # 3. Dump without KMC sorting, then strictly sort using Unix LC_ALL=C sort
        kmc_tools -t{threads} transform "{params.valid_db}" dump "{output.valid_txt}.unsorted"
        LC_ALL=C sort -k1,1 "{output.valid_txt}.unsorted" > "{output.valid_txt}"
        
        # Clean up
        rm -f {params.temp_db}* {params.valid_db}* "{output.valid_txt}.unsorted"
        """

# ------------------------------------------------------------------------------
# Phase 3: Subfamily Filtering & 1% Math
# ------------------------------------------------------------------------------

rule step_07_process_subfamily:
    """Instantly matches local counts against the global master list."""
    input:
        subfam_fasta = os.path.join(config["indexing"]["outputs"]["step_02_fasta_type"], "{repeat}.fasta"),
        valid_txt = config["indexing"]["outputs"]["all_rep_db"] + "_valid.txt"
    output:
        final_fasta = os.path.join(config["indexing"]["outputs"]["step_07_final_fasta"], "{repeat}.fasta")
    params:
        k = config["indexing"]["kmer_len"],
        subfam_db = lambda w, resources: f"{resources.tmpdir}/{w.repeat}_raw",
        mem_gb = lambda wildcards, resources: int(resources.mem_mb / 1000)
    threads: 1 
    resources:
        tmpdir = lambda w: os.path.join(config["paths"]["data"], f"tmp_{w.repeat}"),
        mem_mb = 2000
    container: config["docker_images"]["kmc"]
    shell:
        """
        # 1. Count k-mers in this specific subfamily
        kmc -k{params.k} -m{params.mem_gb} -ci1 -fm -t{threads} "{input.subfam_fasta}" "{params.subfam_db}" {resources.tmpdir}

        # 2. Dump without KMC sorting, then strictly sort using Unix LC_ALL=C sort
        kmc_tools transform "{params.subfam_db}" dump "{params.subfam_db}.unsorted.txt"
        LC_ALL=C sort -k1,1 "{params.subfam_db}.unsorted.txt" > "{params.subfam_db}.txt"

        # 3. Join local and global counts perfectly (Removed the buggy tab flag)
        LC_ALL=C join -1 1 -2 1 "{params.subfam_db}.txt" "{input.valid_txt}" | \
        awk -v name="{wildcards.repeat}" '{{
            if ($2 / $3 >= 0.99) {{
                print ">" name "_" NR "\\n" $1
            }}
        }}' > "{output.final_fasta}"

        # Cleanup local intermediate files
        rm -f {params.subfam_db}* """

# ------------------------------------------------------------------------------
# Phase 4: Reporting
# ------------------------------------------------------------------------------

rule step_08_generate_summary:
    """Generates the final summary of k-mers retained."""
    input:
        files_after = lambda w: expand(
            os.path.join(config["indexing"]["outputs"]["step_07_final_fasta"], "{repeat}.fasta"),
            repeat=get_all_repeat_names(w)
        )
    output:
        summary = config["indexing"]["outputs"]["step_08_summary"]
    params:
        dir_after = config["indexing"]["outputs"]["step_07_final_fasta"]
    shell:
        """
        printf "Name\\tFinal_Kmers\\n" > "{output.summary}"
        for file_path in "{params.dir_after}"/*; do
            filename=$(basename "$file_path")
            count=$(grep -c "^>" "$file_path" || true)
            name="${{filename%.*}}"
            printf "%s\\t%s\\n" "$name" "$count" >> "{output.summary}"
        done
        
        # Cleanup empty FASTAs
        find {params.dir_after} -name "*.fasta" -size 0 -delete
        """
