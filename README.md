# RepVote: Repeat & Gene Quantification Pipeline



RepVote is a specialized Snakemake workflow designed for the concurrent quantification of repetitive elements (via custom Rust tool `repeatquant` / `repeat_counter`) and gene expression (via Kallisto/tximport).

It features a hybrid input system that handles both local FASTQ files and automatic downloading of GSM samples from NCBI SRA, making it easy to integrate public datasets with in-house data.

![RepVote Workflow Diagram](./workflow.svg)

## 🛠️ Standalone `repeatquant` (`repeat_counter`) CLI Usage

The `repeatquant` binary (`repeat_counter`) provides two primary subcommands: `index` and `quant`.

### 1. Indexing (`repeat_counter index`)
Generates filtered subfamily k-mer FASTA files and/or zero-copy binary index files (`.idx`) from genomic reference data.

* **`--genome <PATH>`**: Path to reference genome FASTA file.
* **`--bed <PATH>`**: Path to repeat annotation BED file (0-based: Chr, Start, End, Subfamily#Class/Family).
* **`--gtf <PATH>`**: Path to transcriptome GTF file used to filter out k-mers matching canonical gene transcripts.
* **`--no-transcript-filter`**: Skip transcriptome k-mer filtering.
* **`--from-fasta-dir <PATH>`**: Build zero-copy index directly from pre-filtered FASTA files (conflicts with `--genome`/`--bed`/`--gtf`).
* **`--kmer-fasta-out <PATH>`**: Output directory for inspectable filtered k-mer FASTA files.
* **`--fasta-only`**: Generate inspectable FASTA files only and stop without building binary index.
* **`--output <PATH>`**: Path for output zero-copy index file (`.idx`).
* **`-k, --k <INT>`**: K-mer size (default: 24, supported up to 63).
* **`--threads <INT>`**: Number of parallel threads (default: 8).

### 2. Quantification (`repeat_counter quant`)
Quantifies repeat expression from single-end or paired-end FASTQ files using a compiled zero-copy index.

* **`--index <PATH>`**: Path to zero-copy index (`.idx`) file.
* **`--r1 <PATH>`**: Path to Read 1 FASTQ file (supports gzipped inputs).
* **`--r2 <PATH>`**: Path to Read 2 FASTQ file (optional for paired-end data).
* **`--output <PATH>`**: Path for output TSV count file.
* **`--threads <INT>`**: Number of processing threads (default: 1).
* **`--min-matches <INT>`**: Minimum k-mer matches required to assign a read (default: 1).
* **`--assignments-out <PATH>`**: Optional path to export 1:1 read assignment IDs matching FASTQ order.
* **`--debug`**: Enable debug logging.

---

## 📂 Repository Structure

```text
RepVote
├── workflow/            # Core Snakemake logic (Snakefile, rules, scripts)
├── config/              # Configuration files (samples.tsv, config.yaml)
├── resources/           # Static resources (GTF, BED)
├── data/                # Sequence data (fastq.gz) and results (tsv)
└── benchmarks/          # Performance tracking logs
```

## 🚀 Getting Started

To run RepVote, you can use Docker (easiest, zero-dependency) or a local Conda environment.

### Prerequisites

RepVote requires Snakemake v7.32.4 (newer versions may work but are currently untested).

RepVote use multiple container images via snakemake.

### Installation

Clone the repository:
```text
git clone https://github.com/H-OKW-K/RepVote.git
cd RepVote 
```

## 🧪 Verification / Installation Check

Before running the pipeline on your production data, you can verify your installation and environment setup using the provided dummy files.

1. **Move the dummy resources:** Copy or move the test reference files into the core resources directory:
   
   ```bash
   mv .test/resources/* resources/
   ```
   
2. Configure the test samples: Open config/samples.tsv and add the paths to the dummy FASTQ files located in the .test/fastq/ directory. Also update paths in config/config.yaml as the Configuration section below.

3. Run the pipeline: Execute the workflow using the execution commands outlined below.

## ⚙️ Configuration

Before running the pipeline, you must configure your input samples and file paths in the config/ directory.

### 1. Samples Sheet (config/samples.tsv)

This file defines your input data. The pipeline supports two modes via the type column:

**gsm**: Automatically downloads and processes data from NCBI SRA.

**local**: Uses FASTQ files stored on your disk.

| Column | Description |
| :--- | :--- |
| sample_name | Unique identifier for the sample (Directory name). |
| type | gsm or local. |
| r1 | For GSM: The GSM Accession ID (e.g., GSM123456).<br>For Local: Absolute path to Read 1. |
| r2 | For GSM: Leave empty.<br>For Local: Absolute path to Read 2 (or empty for Single End). |

Example config/samples.tsv:

| sample_name   | type  | r1                          | r2                          |
|:--------------|:------|:----------------------------|:----------------------------|
| Treated_Rep1  | gsm   | GSM1234567                  |                             |
| Control_Rep1  | local | /data/reads/Ctrl_R1.fastq.gz| /data/reads/Ctrl_R2.fastq.gz|
| Control_Rep2  | local | /data/reads/SE_R1.fastq.gz  |                             |

### 2. Pipeline Config (config/config.yaml)
Edit config/config.yaml to set paths for your reference genome and resource limits.

```yaml

---------------------------------------------------------------------

User provided files (REQUIRED)

---------------------------------------------------------------------

resources:
  # The Reference Genome (FASTA)
  genome_fasta: "resources/genome_mock.fa"
  
  # The Repeat Annotation (BED)
  # Format: No Header, 0-based: [Chr, start, end, subfamily#Class/Family]
  repeat_bed: "resources/repeats_mock.bed"
  
  # The Transcriptome (GTF)
  # Used to filter out k-mers that match canonical genes (preventing false positives)
  genome_gtf: "resources/transcripts_mock.gtf"
  
---------------------------------------------------------------------

resource_profiles:
  light:
    threads: 1
    mem_mb: 4000

  standard:
    threads: 4
    mem_mb: 4000
  
  heavy:
    threads: 16      # High threads = Fewer parallel jobs = Safer RAM
    mem_mb: 32000
    
  extreme:
    threads: 24
    mem_mb: 64000

---------------------------------------------------------------------
```

## 🏃 Execution

Ensure Snakemake v7.32.4 or newer is activated in your environment before executing.
```text
    snakemake --cores 16 --use-singularity \
    --resources mem_mb=100000 io_bound_for_esearch=3 io_bound_for_prefetch=3 \
    --rerun-incomplete
```
Note: Customize your runtime resources (e.g., --resources mem_mb=100000) depending on your infrastructure requirements.

## 🧬 Reference Data

Original paper used this pipeline utilized the **T2T-CHM13 v2.0** assembly and corresponding annotations. You can download the compatible versions from the sources below:

| File Type | Description | Source / Download Link |
| :--- | :--- | :--- |
| **Genome** | T2T-CHM13 v2.0 Assembly (`.fa`) | [T2T GitHub (marbl/CHM13)](https://github.com/marbl/CHM13) <br> *File: `chm13v2.0.fa.gz`* |
| **Transcriptome** | NCBI RefSeq Gene Annotations (`.gtf`) | [UCSC GoldenPath (hs1)](https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/genes/) <br> *File: `hs1.ncbiRefSeq.gtf.gz`* |
| **Repeats** | RepeatMasker Annotation (`.bed`) | [T2T GitHub (marbl/CHM13)](https://github.com/marbl/CHM13) <br> *File: `chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.bed`* |
| **Repeats** | Cen/Sat v2.1 Annotation (`.bed`) | [T2T GitHub (marbl/CHM13)](https://github.com/marbl/CHM13) <br> *File: `chm13v2.0_censat_v2.1.bed`* |

## 📊 Outputs

Upon successful completion, results are stored in the results/ directory:

* data/kmer_counts_raw.tsv: Raw count matrix for repetitive elements.
* data/kmer_counts_cpm.tsv: CPM-normalized counts for repeats.
* data/gene_counts_raw.tsv: Gene expression counts (aggregated via tximport).
* data/sample_total_reads.tsv: The number of reads which passed fastp in each library. These value are used for CPM calculation.
* data/MultiQC/multiqc_report.html**: Aggregated quality control report.

## 🤝 Acknowledgements & Credits

### Scientific Inspiration & Methodology
This workflow incorporates a conceptual framework inspired by the **ARTEMIS** pipeline. The methodology of extracting unique k-mers from repeat regions and filtering out transcript-matching k-mers is derived from their work. No scripts from ARTEMIS are directly used in RepVote; all processing is natively implemented in Rust.

* **Paper:** [Annapragada et al., *Science Translational Medicine* (2024)](https://www.science.org/doi/10.1126/scitranslmed.adj9283)
* **Repository:** [cancer-genomics/artemis2024](https://github.com/cancer-genomics/artemis2024)

### Core Tools & Technologies
RepVote stands on the shoulders of giants. We gratefully acknowledge the authors of the following tools integrated into this pipeline:

* **Workflow Management:** [Snakemake](https://snakemake.github.io)
* **Quantification:** Kallisto (Bray et al.) and Tximport (Soneson et al.)
* **QC & Reporting:** FastQC (Andrews), Fastp (Chen et al.), and MultiQC (Ewels et al.)
* **Rust Libraries for repeatquant:** Rayon, Needletail, Rkyv, Memmap2

## 📚 Citation

If you use RepVote in your research, please cite:

> **[Author Name] et al. (Year).** "Title of the Paper." Journal Name, Volume(Issue), Page numbers. DOI: [DOI Link]

## 📜 License

This project is licensed under the MIT License - see the LICENSE file for details.
