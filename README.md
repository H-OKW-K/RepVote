# RepVote: Repeat & Gene Quantification Pipeline



RepVote is a specialized Snakemake workflow designed for the concurrent quantification of repetitive elements (via custom Rust tool `repeatquant` / `repeat_counter`) and gene expression (via Kallisto/tximport).

It features a hybrid input system that handles both local FASTQ files and automatic downloading of GSM samples from NCBI SRA, making it easy to integrate public datasets with in-house data.

![RepVote Workflow Diagram](./workflow.svg)

## 🛠️ Standalone `repeatquant` (`repeat_counter`) CLI Usage

The `repeatquant` binary (`repeat_counter`) provides three primary subcommands: `index`, `quant`, and `characterize`.

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

### 3. Read characterization (`repeat_counter characterize`)
Characterizes reads assigned to selected repeat subfamilies using the assignment outputs produced by `repeat_counter quant`. For paired-end data, R1 and R2 are characterized independently while retaining a shared fragment identifier.

* **`--assignments <PATH>`**: Read-assignment file generated during quantification.
* **`--assignment-map <PATH>`**: Assignment map corresponding to `--assignments`.
* **`--r1 <PATH>`**: Path to Read 1 FASTQ file.
* **`--r2 <PATH>`**: Path to Read 2 FASTQ file (omit for single-end data).
* **`--subfamily <PATH>`**: Path to a text file listing the repeat subfamilies to characterize, one per line. The file must contain the **final repeat names used in the assignment output**, with `/` in the original repeat annotation names replaced by `_` before they are written to the file (e.g., `AluSx#SINE/Alu` becomes `AluSx#SINE_Alu`).
* **`--output-prefix <PATH>`**: Prefix for characterization outputs. The workflow writes `<prefix>.per_read.tsv` and `<prefix>.summary.tsv`.
* **`--threads <INT>`**: Number of processing threads.

#### Per-read characterization columns

| Column | Brief summary | How calculated |
| --- | --- | --- |
| `fragment_id` | Sequential identifier for the assigned fragment. | Incremented once per FASTQ record pair; R1 and R2 from the same paired-end fragment share the same ID. |
| `read_id` | Original FASTQ read identifier. | Taken directly from the FASTQ header for the individual read. |
| `mate` | Indicates which observed read was characterized. | `R1` for read 1 and `R2` for read 2; each mate is characterized independently. |
| `subfamily` | User-selected subfamily associated with this characterized read. | Intersection of the fragment's winning assignments with the subfamilies requested for characterization; multiple names are comma-separated if ambiguous assignments are included. |
| `assignment` | Full repeat-subfamily assignment of the fragment. | R1 and R2 k-mer votes are combined; all subfamilies tied for the highest vote count above `min_matches` are reported. |
| `assignment_ambiguous` | Indicates whether the fragment has multiple equally supported assignments. | `true` when more than one subfamily ties for the maximum assignment score; otherwise `false`. |
| `length` | Total observed read length in nucleotides. | Number of sequence characters in the FASTQ read, including ambiguous bases such as `N`. |
| `valid_bases` | Number of bases with definite nucleotide identity. | Count of positions containing A, C, G, or T; ambiguous bases are excluded. |
| `gc_fraction` | Fraction of valid bases that are G or C. | `(G + C) / valid_bases`. |
| `a_fraction` | Fraction of valid bases that are A. | `A / valid_bases`. |
| `c_fraction` | Fraction of valid bases that are C. | `C / valid_bases`. |
| `g_fraction` | Fraction of valid bases that are G. | `G / valid_bases`. |
| `t_fraction` | Fraction of valid bases that are T. | `T / valid_bases`. |
| `base_entropy` | Normalized nucleotide-composition entropy; higher values indicate more balanced base composition. | Shannon entropy `−Σ pᵢ log₂(pᵢ)` over A/C/G/T frequencies, divided by 2 so the range is 0–1. |
| `max_homopolymer_nt` | Length of the longest consecutive run of one nucleotide. | Maximum run of identical A, C, G, or T bases; ambiguous bases break a run. |
| `mono_dominance` | Fraction of valid 1-mer positions occupied by the most frequent nucleotide. | Count of the most frequent 1-mer divided by the number of valid 1-mer positions. |
| `di_dominance` | Fraction of valid 2-mer positions occupied by the most frequent dinucleotide. | Count of the most frequent 2-mer divided by all valid overlapping 2-mer windows. Windows containing ambiguous bases are excluded. |
| `tri_dominance` | Fraction of valid 3-mer positions occupied by the most frequent trinucleotide. | Count of the most frequent 3-mer divided by all valid overlapping 3-mer windows. Windows containing ambiguous bases are excluded. |
| `rc_covered_fraction` | Fraction of valid bases participating in qualifying exact intraread reverse-complementary tracts. | Positions belonging to ≥1 exact reverse-complementary tract of at least `min_rc_tract` nt are marked; marked positions are divided by `valid_bases`. |
| `max_rc_tract_nt` | Length of the longest exact reverse-complementary tract within the read. | Finds maximal exact reverse-complementary alignments between distinct regions of the read and reports the longest qualifying tract. |
| `rc_tract_count` | Number of qualifying exact reverse-complementary tracts. | Counts maximal reverse-complementary alignment runs of at least `min_rc_tract` nt, avoiding counting overlapping k-mer matches from the same alignment as separate events. |
| `dominant_period_nt` | Shortest strongly supported tandem-repeat period. | Tests periods `1 … min(max_period, length/2)` by comparing bases separated by each period; among near-best qualifying periods, the shortest period is selected. `NA` if none reaches the threshold. |
| `periodicity_score` | Strength of sequence agreement with the selected repeat period. | For period `p`, number of valid position pairs where `base[i] = base[i−p]` divided by the number of comparable valid pairs. If no period qualifies, the best observed score is reported. |
| `repeat_units` | Approximate number of repeat units represented by the read. | `length / dominant_period_nt`; reported as `NA` when no dominant period is detected. |

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
* data/<sample>/characterize/selected.per_read.tsv: Per-read characterization metrics for the configured repeat subfamilies (when assignment export/characterization is enabled).
* data/<sample>/characterize/selected.summary.tsv: Summary characterization output for the configured repeat subfamilies.
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
