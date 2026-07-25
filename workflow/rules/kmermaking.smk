rule generate_final_kmers_and_index:
    """
    Generates final filtered subfamily FASTA files AND builds the zero-copy binary index.
    """
    container: config["docker_images"]["repeatquant"]
    input:
        genome = config["resources"]["genome_fasta"],
        bed = config["resources"]["repeat_bed"],
        gtf = config["resources"]["genome_gtf"]
    output:
        index = config["quantification"]["repeat_index_path"],
        fasta_dir = directory(config["indexing"]["outputs"]["final_fasta"])
    params:
        k = config["indexing"]["kmer_len"],
        no_transcript_filter = "--no-transcript-filter" if config["indexing"].get("no_transcript_filter", False) else ""
    threads: config["resource_profiles"]["extreme"]["threads"]
    shell:
        """
        repeat_counter index \
            --genome {input.genome} \
            --bed {input.bed} \
            --gtf {input.gtf} \
            --kmer-fasta-out {output.fasta_dir} \
            --output {output.index} \
            --k {params.k} \
            {params.no_transcript_filter} \
            --threads {threads}
        """

rule generate_transcriptome:
    """Generates standard transcriptome FASTA for Kallisto quantification."""
    container: config["docker_images"]["gffread"]
    input:
        genome = config["resources"]["genome_fasta"],
        gtf = config["resources"]["genome_gtf"]
    output:
        fasta = config["indexing"]["generated_transcriptome"]
    shell:
        "gffread -w '{output.fasta}' -g '{input.genome}' '{input.gtf}'"
