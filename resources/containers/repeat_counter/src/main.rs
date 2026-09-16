use clap::{Parser, Subcommand};
use needletail::{parse_fastx_reader, FastxReader};
use rayon::prelude::*;
use rkyv::{Archive, Deserialize, Serialize};
use memmap2::MmapOptions;
use rustc_hash::{FxHashMap, FxHashSet};
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::path::PathBuf;
use std::process::{Command, Stdio};
use std::sync::mpsc::sync_channel;
use std::thread;

#[derive(Parser)]
#[command(name = "repeat_counter", version = "0.3.1")]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Generate final filtered k-mer FASTAs and/or zero-copy binary index
    Index {
        /// Reference genome FASTA
        #[arg(long, requires = "bed")]
        genome: Option<PathBuf>,

        /// Repeat annotation BED file
        #[arg(long)]
        bed: Option<PathBuf>,

        /// Transcriptome GTF file (Optional)
        #[arg(long)]
        gtf: Option<PathBuf>,

        /// Skip transcriptome k-mer filtering
        #[arg(long, default_value_t = false)]
        no_transcript_filter: bool,

        /// Direct input directory containing pre-filtered FASTA files
        #[arg(long, conflicts_with_all = ["genome", "bed", "gtf"])]
        from_fasta_dir: Option<PathBuf>,

        /// Directory to output user-inspectable filtered k-mer FASTA files
        #[arg(long)]
        kmer_fasta_out: Option<PathBuf>,

        /// Stop after generating inspectable FASTA files (Do not build zero-copy index)
        #[arg(long, requires = "kmer_fasta_out")]
        fasta_only: bool,

        /// Output zero-copy index file (.idx)
        #[arg(long, required_unless_present = "fasta_only")]
        output: Option<PathBuf>,

        /// K-mer size (supported up to 63)
        #[arg(short, long, default_value_t = 24)]
        k: usize,

        /// Number of parallel threads
        #[arg(long, default_value_t = 8)]
        threads: usize,
    },

    /// Quantify repeat expression from FASTQ files
    Quant {
        #[arg(long)]
        index: PathBuf,
        #[arg(long)]
        r1: PathBuf,
        #[arg(long)]
        r2: Option<PathBuf>,
        #[arg(long)]
        output: PathBuf,
        #[arg(long, default_value_t = 1)]
        threads: usize,
        /// Minimum k-mer matches required to assign a read (default: 1)
        #[arg(long, default_value_t = 1)]
        min_matches: u32,
        /// Optional path to export 1:1 read assignment IDs matching FASTQ order
        #[arg(long)]
        assignments_out: Option<PathBuf>,
        #[arg(long)]
        debug: bool,
    },

    /// Characterize sequence properties of reads assigned to selected repeat subfamilies
    Characterize {
        /// FASTQ-order fragment assignments exported by quant --assignments-out
        #[arg(long)]
        assignments: PathBuf,
        /// Optional assignment ID-to-subfamily map. Defaults to <assignments>.map.tsv
        #[arg(long)]
        assignment_map: Option<PathBuf>,
        /// FASTQ read 1 (or single-end FASTQ)
        #[arg(long)]
        r1: PathBuf,
        /// Optional paired-end FASTQ read 2
        #[arg(long)]
        r2: Option<PathBuf>,
        /// Repeat subfamily to characterize. May be repeated or comma-delimited.
        #[arg(long, required = true, value_delimiter = ',')]
        subfamily: Vec<String>,
        /// Output prefix. Produces <prefix>.per_read.tsv and <prefix>.summary.tsv
        #[arg(long)]
        output_prefix: PathBuf,
        /// Number of parallel threads
        #[arg(long, default_value_t = 1)]
        threads: usize,
        /// Include tied assignments that contain a selected subfamily
        #[arg(long, default_value_t = false)]
        include_ambiguous: bool,
        /// Minimum exact reverse-complement tract length to retain
        #[arg(long, default_value_t = 5)]
        min_rc_tract: usize,
        /// Maximum tandem-repeat period to test
        #[arg(long, default_value_t = 12)]
        max_period: usize,
        /// Minimum periodicity score required to report a dominant period
        #[arg(long, default_value_t = 0.80)]
        min_periodicity_score: f64,
        /// Prefer a shorter period when its score is within this delta of the best score
        #[arg(long, default_value_t = 0.02)]
        period_equivalence_delta: f64,
        #[arg(long)]
        debug: bool,
    },
}

#[derive(Archive, Deserialize, Serialize)]
#[archive(check_bytes)]
pub struct Partition {
    pub map: FxHashMap<u128, Vec<u16>>,
}

#[derive(Archive, Deserialize, Serialize)]
#[archive(check_bytes)]
pub struct Metadata {
    pub k: u8,
    pub names: Vec<String>,
    pub offsets: Vec<(u64, u64)>,
}

fn log_memory_usage(stage: &str) {
    if let Ok(statm) = fs::read_to_string("/proc/self/statm") {
        let parts: Vec<&str> = statm.split_whitespace().collect();
        if parts.len() >= 2 {
            if let (Ok(virt_pages), Ok(rss_pages)) =
                (parts[0].parse::<usize>(), parts[1].parse::<usize>())
            {
                let page_size_kb = 4; // standard 4KB page size
                let rss_mb = (rss_pages * page_size_kb) / 1024;
                let virt_mb = (virt_pages * page_size_kb) / 1024;
                eprintln!(
                    "[MEM LOG] Stage: '{}' | RSS: {} MB | VIRT: {} MB",
                    stage, rss_mb, virt_mb
                );
            }
        }
    }
}

fn main() {
    let args = Cli::parse();

    match args.command {
        Commands::Index {
            genome,
            bed,
            gtf,
            no_transcript_filter,
            from_fasta_dir,
            kmer_fasta_out,
            fasta_only,
            output,
            k,
            threads,
        } => {
            rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build_global()
                .unwrap();

            if let Some(fasta_dir) = from_fasta_dir {
                let out_idx =
                    output.expect("--output is required when building from fasta directory");
                build_index_from_fasta_dir(&fasta_dir, &out_idx, k as u8);
            } else {
                build_full_pipeline(
                    genome.expect("--genome is required when building from BED"),
                    bed.expect("--bed is required when building from BED"),
                    gtf,
                    no_transcript_filter,
                    kmer_fasta_out,
                    fasta_only,
                    output,
                    k as u8,
                );
            }
        }
        Commands::Quant {
            index,
            r1,
            r2,
            output,
            threads,
            min_matches,
            assignments_out,
            debug,
        } => {
            rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build_global()
                .unwrap();
            quantify_sample(index, r1, r2, output, min_matches, assignments_out, debug);
        }
        Commands::Characterize {
            assignments,
            assignment_map,
            r1,
            r2,
            subfamily,
            output_prefix,
            threads,
            include_ambiguous,
            min_rc_tract,
            max_period,
            min_periodicity_score,
            period_equivalence_delta,
            debug,
        } => {
            rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build_global()
                .unwrap();
            characterize_sample(
                assignments,
                assignment_map,
                r1,
                r2,
                subfamily,
                output_prefix,
                include_ambiguous,
                min_rc_tract,
                max_period,
                min_periodicity_score,
                period_equivalence_delta,
                debug,
            );
        }
    }
}

// ==============================================================================
// INDEXING LOGIC
// ==============================================================================

fn build_index_from_fasta_dir(fasta_dir: &PathBuf, output: &PathBuf, k: u8) {
    log_memory_usage("Start build_index_from_fasta_dir");
    println!(
        "[INFO] Reading final filtered FASTA files from {:?}...",
        fasta_dir
    );

    let mut fasta_files: Vec<PathBuf> = fs::read_dir(fasta_dir)
        .expect("Cannot read FASTA directory")
        .filter_map(|e| e.ok().map(|x| x.path()))
        .filter(|p| p.extension().map_or(false, |ext| ext == "fasta" || ext == "fa"))
        .collect();

    fasta_files.sort();

    let mut temp_dir = output.clone();
    temp_dir.set_extension("tmp_partitions");
    fs::create_dir_all(&temp_dir).expect("Failed to create temp directory");

    let mut writers: Vec<BufWriter<File>> = (0..256)
        .map(|i| {
            let p = temp_dir.join(format!("part_{:03}.bin", i));
            BufWriter::with_capacity(64 * 1024, File::create(p).unwrap())
        })
        .collect();

    let mut names: Vec<String> = Vec::new();

    for (family_id, path) in fasta_files.iter().enumerate() {
        let stem = path.file_stem().unwrap().to_str().unwrap().to_string();
        names.push(stem);

        let mut reader = needletail::parse_fastx_file(path).expect("Invalid FASTA file");
        while let Some(Ok(rec)) = reader.next() {
            let seq = rec.seq();
            if seq.len() == k as usize {
                if let Some(enc) = encode_kmer(&seq) {
                    let rc_enc = get_rc_kmer(enc, k as usize);
                    let canonical = std::cmp::min(enc, rc_enc);
                    let part_idx = (canonical & 0xFF) as usize;

                    writers[part_idx]
                        .write_all(&canonical.to_le_bytes())
                        .unwrap();
                    writers[part_idx]
                        .write_all(&(family_id as u16).to_le_bytes())
                        .unwrap();
                }
            }
        }
    }

    for w in writers.iter_mut() {
        w.flush().ok();
    }
    drop(writers);

    log_memory_usage("Before write_zero_copy_index from fasta_dir");
    println!("[INFO] Assembling zero-copy index...");
    write_zero_copy_index(output, &temp_dir, names, k);
    fs::remove_dir_all(temp_dir).ok();
    log_memory_usage("Complete build_index_from_fasta_dir");
    println!("[INFO] Indexing complete.");
}

struct BedRecord {
    start: usize,
    end: usize,
    subfam_id: usize,
}

fn get_family_from_name(name: &str) -> &str {
    let subfam = if let Some((_, after)) = name.split_once('#') {
        after
    } else {
        name
    };
    if let Some((before, _)) = subfam.split_once('/') {
        before
    } else {
        subfam
    }
}

fn sanitize_name(name: &str) -> String {
    let s = name.replace('/', "_");
    s.chars()
        .filter(|c| !matches!(c, '?' | '*' | '<' | '>' | '|' | ':'))
        .collect()
}

fn decode_kmer(val: u128, k: u8) -> String {
    let mut bytes = vec![0u8; k as usize];
    for i in 0..k as usize {
        let shift = 2 * (k as usize - 1 - i);
        let base = (val >> shift) & 3;
        bytes[i] = match base {
            0 => b'A',
            1 => b'C',
            2 => b'G',
            3 => b'T',
            _ => unreachable!(),
        };
    }
    String::from_utf8(bytes).unwrap()
}

fn build_full_pipeline(
    genome_path: PathBuf,
    bed_path: PathBuf,
    gtf_path: Option<PathBuf>,
    no_transcript_filter: bool,
    kmer_fasta_out: Option<PathBuf>,
    fasta_only: bool,
    output_index: Option<PathBuf>,
    k: u8,
) {
    log_memory_usage("Start build_full_pipeline");
    println!("[INFO] Running full filter pipeline with disk spilling & u128 k-mers (Genome, Repeat, Transcript subtraction)...");

    let fasta_dir = kmer_fasta_out.unwrap_or_else(|| PathBuf::from("results/kmer_fastas"));
    fs::create_dir_all(&fasta_dir).expect("Failed to create kmer fasta output dir");

    let temp_dir = fasta_dir.join("tmp_partitions");
    fs::create_dir_all(&temp_dir).expect("Failed to create temporary partitions directory");

    // 1. Parse BED file
    println!("[INFO] Parsing repeat BED file {:?}...", bed_path);
    let bed_file = File::open(&bed_path).expect("Failed to open BED file");
    let bed_reader = BufReader::new(bed_file);

    let mut subfam_names: Vec<String> = Vec::new();
    let mut subfam_map: FxHashMap<String, usize> = FxHashMap::default();
    let mut bed_records_by_chrom: FxHashMap<String, Vec<BedRecord>> = FxHashMap::default();
    let mut repeat_intervals_by_chrom: FxHashMap<String, Vec<(usize, usize)>> = FxHashMap::default();

    for line_res in bed_reader.lines() {
        let line = match line_res {
            Ok(l) => l,
            Err(_) => continue,
        };
        let trimmed = line.trim();
        if trimmed.is_empty()
            || trimmed.starts_with('#')
            || trimmed.starts_with("track")
            || trimmed.starts_with("browser")
        {
            continue;
        }
        let cols: Vec<&str> = trimmed.split('\t').collect();
        if cols.len() < 4 {
            continue;
        }

        let chrom = cols[0].to_string();
        let start: usize = match cols[1].parse() {
            Ok(s) => s,
            Err(_) => continue,
        };
        let end: usize = match cols[2].parse() {
            Ok(e) => e,
            Err(_) => continue,
        };

        if end <= start {
            continue;
        }

        let raw_name = cols[3];
        let family = get_family_from_name(raw_name);
        if family == "Unknown" || family == "Simple_repeat" || family == "Low_complexity" {
            continue;
        }
       
        repeat_intervals_by_chrom
            .entry(chrom.clone())
            .or_default()
            .push((start, end));

        let san_name = sanitize_name(raw_name);
        let subfam_id = *subfam_map.entry(san_name.clone()).or_insert_with(|| {
            let id = subfam_names.len();
            subfam_names.push(san_name);
            id
        });

        bed_records_by_chrom
            .entry(chrom)
            .or_default()
            .push(BedRecord {
                start,
                end,
                subfam_id,
            });
    }
    log_memory_usage("After BED parsing");

    // 2. Parse GTF file (if transcript filtering enabled)
    struct Exon {
        start: usize,
        end: usize,
    }
    let mut tx_exons_by_chrom: FxHashMap<String, Vec<(String, Exon)>> = FxHashMap::default();
    let use_gtf = !no_transcript_filter && gtf_path.as_ref().map_or(false, |p| p.exists());

    if use_gtf {
        let gtf_p = gtf_path.as_ref().unwrap();
        println!("[INFO] Parsing GTF annotations {:?}...", gtf_p);
        let gtf_file = File::open(gtf_p).expect("Cannot open GTF file");
        let reader = BufReader::new(gtf_file);

        for line_res in reader.lines() {
            let line = match line_res {
                Ok(l) => l,
                Err(_) => continue,
            };
            if line.starts_with('#') {
                continue;
            }
            let cols: Vec<&str> = line.split('\t').collect();
            if cols.len() < 9 || cols[2] != "exon" {
                continue;
            }
            let chrom = cols[0].to_string();
            let start: usize = match cols[3].parse::<usize>() {
                Ok(s) if s > 0 => s - 1,
                _ => continue,
            };
            let end: usize = match cols[4].parse() {
                Ok(e) => e,
                _ => continue,
            };

            let attr = cols[8];
            let tx_id = attr.split(';').find_map(|kv| {
                let kv = kv.trim();
                if kv.starts_with("transcript_id") {
                    let parts: Vec<&str> = kv.split_whitespace().collect();
                    if parts.len() >= 2 {
                        Some(parts[1].trim_matches('"').to_string())
                    } else {
                        None
                    }
                } else {
                    None
                }
            });

            if let Some(id) = tx_id {
                tx_exons_by_chrom
                    .entry(chrom)
                    .or_default()
                    .push((id, Exon { start, end }));
            }
        }
        log_memory_usage("After GTF parsing");
    }

    // Initialize 64KB buffered partition writers for chunked flushing
    let mut kg_writers: Vec<BufWriter<File>> = (0..256)
        .map(|i| {
            let p = temp_dir.join(format!("kg_{:03}.bin", i));
            BufWriter::with_capacity(64 * 1024, File::create(p).unwrap())
        })
        .collect();

    let mut repeat_writers: Vec<BufWriter<File>> = (0..256)
        .map(|i| {
            let p = temp_dir.join(format!("repeat_{:03}.bin", i));
            BufWriter::with_capacity(64 * 1024, File::create(p).unwrap())
        })
        .collect();

    let mut kt_writers: Vec<BufWriter<File>> = (0..256)
        .map(|i| {
            let p = temp_dir.join(format!("kt_{:03}.bin", i));
            BufWriter::with_capacity(64 * 1024, File::create(p).unwrap())
        })
        .collect();

    let k_usize = k as usize;
    let mut tx_builder_map: FxHashMap<String, Vec<(usize, Vec<u8>)>> = FxHashMap::default();

    // 3. Stream Reference Genome FASTA chromosome by chromosome
    println!("[INFO] Streaming genome FASTA and spilling KG & Repeat k-mers to 256 disk partitions...");
    let mut reader = needletail::parse_fastx_file(&genome_path).expect("Failed to open genome FASTA");

    let mut chrom_count = 0;
    while let Some(Ok(rec)) = reader.next() {
        let header = std::str::from_utf8(rec.id()).unwrap();
        let chrom = header.split_whitespace().next().unwrap().to_string();
        let seq = rec.seq().to_ascii_uppercase();
        let len = seq.len();
        chrom_count += 1;

        // 3a. Mask unmasked genome k-mers (Set KG)
        let mut mask = vec![false; len];
        if let Some(intervals) = repeat_intervals_by_chrom.get(&chrom) {
            for &(start, end) in intervals {
                let e = std::cmp::min(end, len);
                if start < e {
                    for p in start..e {
                        mask[p] = true;
                    }
                }
            }
        }

        if len >= k_usize {
            let mut run = 0;
            for i in 0..len {
                if mask[i] {
                    run = 0;
                } else {
                    run += 1;
                    if run >= k_usize {
                        let start_idx = i + 1 - k_usize;
                        if let Some(enc) = encode_kmer(&seq[start_idx..=i]) {
                            let rc = get_rc_kmer(enc, k_usize);
                            let canonical = std::cmp::min(enc, rc);
                            let part_idx = (canonical & 0xFF) as usize;
                            kg_writers[part_idx]
                                .write_all(&canonical.to_le_bytes())
                                .unwrap();
                        }
                    }
                }
            }
        }

        // 3b. Stream repeat k-mer instances with subfamily ID
        if let Some(records) = bed_records_by_chrom.get(&chrom) {
            for record in records {
                let e = std::cmp::min(record.end, len);
                if record.start < e && e - record.start >= k_usize {
                    let sub_seq = &seq[record.start..e];
                    for i in 0..=(sub_seq.len() - k_usize) {
                        if let Some(enc) = encode_kmer(&sub_seq[i..i + k_usize]) {
                            let rc = get_rc_kmer(enc, k_usize);
                            let canonical = std::cmp::min(enc, rc);
                            let part_idx = (canonical & 0xFF) as usize;
                            repeat_writers[part_idx]
                                .write_all(&canonical.to_le_bytes())
                                .unwrap();
                            repeat_writers[part_idx]
                                .write_all(&(record.subfam_id as u16).to_le_bytes())
                                .unwrap();
                        }
                    }
                }
            }
        }

        // 3c. Collect GTF exons for transcriptome k-mers (Set KT)
        if use_gtf {
            if let Some(exons) = tx_exons_by_chrom.get(&chrom) {
                for (tx_id, exon) in exons {
                    let e = std::cmp::min(exon.end, len);
                    if exon.start < e {
                        let slice = seq[exon.start..e].to_vec();
                        tx_builder_map
                            .entry(tx_id.clone())
                            .or_default()
                            .push((exon.start, slice));
                    }
                }
            }
        }

        if chrom_count % 5 == 0 {
            log_memory_usage(&format!("Streaming genome (chrom {})", chrom_count));
        }
    }

    // Flush KG and Repeat partition files
    for w in kg_writers.iter_mut() {
        w.flush().ok();
    }
    for w in repeat_writers.iter_mut() {
        w.flush().ok();
    }
    log_memory_usage("After streaming genome & flushing KG/Repeat partitions");

    // 4. Build and stream transcriptome k-mers (Set KT) in chunks
    if use_gtf {
        println!("[INFO] Processing spliced transcript k-mers (Set KT)...");
        for (_tx_id, mut parts) in tx_builder_map {
            parts.sort_by_key(|p| p.0);
            let mut spliced = Vec::new();
            for (_, slice) in parts {
                spliced.extend_from_slice(&slice);
            }
            if spliced.len() >= k_usize {
                for i in 0..=(spliced.len() - k_usize) {
                    if let Some(enc) = encode_kmer(&spliced[i..i + k_usize]) {
                        let rc = get_rc_kmer(enc, k_usize);
                        let canonical = std::cmp::min(enc, rc);
                        let part_idx = (canonical & 0xFF) as usize;
                        kt_writers[part_idx]
                            .write_all(&canonical.to_le_bytes())
                            .unwrap();
                    }
                }
            }
        }
    }

    for w in kt_writers.iter_mut() {
        w.flush().ok();
    }
    drop(kg_writers);
    drop(repeat_writers);
    drop(kt_writers);
    log_memory_usage("After flushing KT partitions");

    // 5. Parallel Partition Processing & Cache-Friendly RLE
    println!("[INFO] Processing 256 partitions in parallel (Cache-friendly RLE & fast set filtering)...");
    let num_subfams = subfam_names.len();

    (0..256)
        .into_par_iter()
        .for_each(|part_idx| {
            let kg_path = temp_dir.join(format!("kg_{:03}.bin", part_idx));
            let kt_path = temp_dir.join(format!("kt_{:03}.bin", part_idx));
            let rep_path = temp_dir.join(format!("repeat_{:03}.bin", part_idx));

            // Load and sort KG
            let mut kg_sorted: Vec<u128> = Vec::new();
            if kg_path.exists() {
                if let Ok(mut f) = File::open(&kg_path) {
                    let mut buf = [0u8; 16];
                    while f.read_exact(&mut buf).is_ok() {
                        kg_sorted.push(u128::from_le_bytes(buf));
                    }
                }
                kg_sorted.sort_unstable();
                kg_sorted.dedup();
            }

            // Load and sort KT
            let mut kt_sorted: Vec<u128> = Vec::new();
            if kt_path.exists() {
                if let Ok(mut f) = File::open(&kt_path) {
                    let mut buf = [0u8; 16];
                    while f.read_exact(&mut buf).is_ok() {
                        kt_sorted.push(u128::from_le_bytes(buf));
                    }
                }
                kt_sorted.sort_unstable();
                kt_sorted.dedup();
            }

            // Load and sort Repeat entries
            let mut repeat_entries: Vec<(u128, u16)> = Vec::new();
            if rep_path.exists() {
                if let Ok(mut f) = File::open(&rep_path) {
                    let mut buf = [0u8; 18];
                    while f.read_exact(&mut buf).is_ok() {
                        let kmer = u128::from_le_bytes(buf[0..16].try_into().unwrap());
                        let subfam = u16::from_le_bytes(buf[16..18].try_into().unwrap());
                        repeat_entries.push((kmer, subfam));
                    }
                }
                repeat_entries.sort_unstable_by_key(|&(kmer, _)| kmer);
            }

            let part_out_path = temp_dir.join(format!("part_{:03}.bin", part_idx));
            let mut part_writer =
                BufWriter::with_capacity(64 * 1024, File::create(part_out_path).unwrap());

            let mut i = 0;
            while i < repeat_entries.len() {
                let current_kmer = repeat_entries[i].0;
                let start = i;

                while i < repeat_entries.len() && repeat_entries[i].0 == current_kmer {
                    i += 1;
                }
                let c_global = (i - start) as u32;

                // Exact Set Check via O(log N) Binary Search
                if kg_sorted.binary_search(&current_kmer).is_ok()
                    || kt_sorted.binary_search(&current_kmer).is_ok()
                {
                    continue;
                }

                // Count C_local per subfamily in contiguous slice
                let mut subfam_counts: FxHashMap<u16, u32> = FxHashMap::default();
                for entry in &repeat_entries[start..i] {
                    *subfam_counts.entry(entry.1).or_insert(0) += 1;
                }

                for (subfam_id, c_local) in subfam_counts {
                    if (c_local as f64) / (c_global as f64) >= 0.99 {
                        part_writer.write_all(&current_kmer.to_le_bytes()).unwrap();
                        part_writer.write_all(&(subfam_id as u16).to_le_bytes()).unwrap();
                    }
                }
            }
            part_writer.flush().ok();
        });

    log_memory_usage("After 256 partition RLE filtering");

    // 6. Write Filtered FASTAs
    let mut fasta_writers: Vec<Option<BufWriter<File>>> = (0..num_subfams).map(|_| None).collect();
    let mut subfam_kmer_counts = vec![0usize; num_subfams];

    for part_idx in 0..256 {
        let p = temp_dir.join(format!("part_{:03}.bin", part_idx));
        if p.exists() {
            if let Ok(mut f) = File::open(&p) {
                let mut buf = [0u8; 18];
                while f.read_exact(&mut buf).is_ok() {
                    let kmer = u128::from_le_bytes(buf[0..16].try_into().unwrap());
                    let subfam_id = u16::from_le_bytes(buf[16..18].try_into().unwrap()) as usize;
                    subfam_kmer_counts[subfam_id] += 1;
                    let count = subfam_kmer_counts[subfam_id];
                    let san_name = &subfam_names[subfam_id];
                    let writer = fasta_writers[subfam_id].get_or_insert_with(|| {
                        let out_path = fasta_dir.join(format!("{}.fasta", san_name));
                        BufWriter::new(File::create(out_path).expect("Failed to create fasta file"))
                    });
                    let seq_str = decode_kmer(kmer, k);
                    writeln!(writer, ">{}_{}\n{}", san_name, count, seq_str).unwrap();
                }
            }
        }
    }
    for w in fasta_writers.iter_mut().flatten() {
        w.flush().ok();
    }

    println!("[INFO] Final filtered FASTAs written to {:?}", fasta_dir);

    if !fasta_only {
        let out_idx = output_index
            .expect("--output path is required when --fasta-only is not present");
        log_memory_usage("Before zero-copy index assembly");

        write_zero_copy_index(&out_idx, &temp_dir, subfam_names, k);
    }

    fs::remove_dir_all(temp_dir).ok();
    log_memory_usage("Complete build_full_pipeline");
}

fn write_zero_copy_index(output: &PathBuf, temp_dir: &PathBuf, names: Vec<String>, k: u8) {
    let mut out_writer = BufWriter::new(File::create(output).unwrap());
    let mut offsets = Vec::with_capacity(256);
    let mut current_offset: u64 = 0;

    for i in 0..256 {
        let p = temp_dir.join(format!("part_{:03}.bin", i));
        let mut map: FxHashMap<u128, Vec<u16>> = FxHashMap::default();

        if p.exists() {
            let mut f = File::open(&p).unwrap();
            let file_len = f.metadata().map(|m| m.len()).unwrap_or(0);
            let estimated_entries = (file_len / 18) as usize;
            map.reserve(estimated_entries);

            let mut buf = [0u8; 18];
            while f.read_exact(&mut buf).is_ok() {
                let kmer = u128::from_le_bytes(buf[0..16].try_into().unwrap());
                let id = u16::from_le_bytes(buf[16..18].try_into().unwrap());
                map.entry(kmer).or_default().push(id);
            }
            fs::remove_file(p).ok();
        }

        let partition = Partition { map };
        let bytes = rkyv::to_bytes::<_, 4096>(&partition).unwrap();
        out_writer.write_all(&bytes).unwrap();

        let chunk_len = bytes.len() as u64;
        offsets.push((current_offset, current_offset + chunk_len));
        current_offset += chunk_len;
    }

    let meta = Metadata { k, names, offsets };
    let meta_bytes = rkyv::to_bytes::<_, 4096>(&meta).unwrap();
    out_writer.write_all(&meta_bytes).unwrap();
    out_writer.write_all(&current_offset.to_le_bytes()).unwrap();
}

// ==============================================================================
// QUANTIFICATION & READ ASSIGNMENT LOGIC
// ==============================================================================

fn quantify_sample(
    index_path: PathBuf,
    r1_path: PathBuf,
    r2_path: Option<PathBuf>,
    output_path: PathBuf,
    min_matches: u32,
    assignments_out: Option<PathBuf>,
    debug: bool,
) {
    if debug {
        eprintln!(
            "[INFO] Loading zero-copy index with min_matches = {}...",
            min_matches
        );
    }

    let file = File::open(&index_path).expect("Index not found");
    let mmap = unsafe { MmapOptions::new().map(&file).expect("Failed to mmap") };
    let total_len = mmap.len();

    let meta_ptr_bytes: [u8; 8] = mmap[total_len - 8..total_len].try_into().unwrap();
    let meta_offset = u64::from_le_bytes(meta_ptr_bytes) as usize;

    let meta = unsafe { rkyv::archived_root::<Metadata>(&mmap[meta_offset..total_len - 8]) };
    let k = meta.k as usize;
    let num_classes = meta.names.len();

    let mut partitions = Vec::with_capacity(256);
    for i in 0..256 {
        let start = meta.offsets[i].0 as usize;
        let end = meta.offsets[i].1 as usize;
        let archived_part = unsafe { rkyv::archived_root::<Partition>(&mmap[start..end]) };
        partitions.push(&archived_part.map);
    }

    // Optional Read Assignment Output Writers
    let mut assign_writer = assignments_out.as_ref().map(|path| {
        let map_path = default_assignment_map_path(path);

        let mut map_file = BufWriter::new(File::create(map_path).expect("Cannot create ID map file"));
        writeln!(map_file, "ID\tSubfamily").unwrap();
        for (id, name) in meta.names.iter().enumerate() {
            writeln!(map_file, "{}\t{}", id, name).unwrap();
        }

        BufWriter::new(File::create(path).expect("Cannot create assignments output file"))
    });

    let mut final_counts = vec![0.0; num_classes];
    let (tx, rx) = sync_channel::<Vec<(Vec<u8>, Option<Vec<u8>>)>>(4);

    thread::spawn(move || {
        let batch_size = 50_000;
        let mut r1 = get_fastx_reader(&r1_path);
        if let Some(r2_p) = r2_path {
            let mut r2 = get_fastx_reader(&r2_p);
            loop {
                let mut batch = Vec::with_capacity(batch_size);
                for _ in 0..batch_size {
                    match (r1.next(), r2.next()) {
                        (Some(Ok(a)), Some(Ok(b))) => {
                            batch.push((a.seq().to_vec(), Some(b.seq().to_vec())))
                        }
                        _ => break,
                    }
                }
                if batch.is_empty() {
                    break;
                }
                if tx.send(batch).is_err() {
                    break;
                }
            }
        } else {
            loop {
                let mut batch = Vec::with_capacity(batch_size);
                for _ in 0..batch_size {
                    if let Some(Ok(a)) = r1.next() {
                        batch.push((a.seq().to_vec(), None));
                    } else {
                        break;
                    }
                }
                if batch.is_empty() {
                    break;
                }
                if tx.send(batch).is_err() {
                    break;
                }
            }
        }
    });

    for batch in rx {
        let (batch_counts, batch_assignments) = process_batch_rkyv(
            &batch,
            k,
            &partitions,
            num_classes,
            min_matches,
            assignments_out.is_some(),
        );

        for (i, c) in batch_counts.iter().enumerate() {
            final_counts[i] += c;
        }

        if let (Some(writer), Some(assignments)) = (&mut assign_writer, batch_assignments) {
            for assign_id in assignments {
                writeln!(writer, "{}", assign_id).unwrap();
            }
        }
    }

    let mut out = File::create(output_path).unwrap();
    writeln!(out, "RepeatType\tCount").unwrap();
    for (i, count) in final_counts.iter().enumerate() {
        writeln!(out, "{}\t{:.2}", meta.names[i], count).unwrap();
    }
}


// ============================================================================== 
// POST-ASSIGNMENT READ CHARACTERIZATION
// ============================================================================== 

#[derive(Clone)]
struct CharacterizationFragment {
    fragment_id: u64,
    r1_id: String,
    r1_seq: Vec<u8>,
    r2: Option<(String, Vec<u8>)>,
    winners: Vec<u16>,
    selected_winners: Vec<u16>,
}

#[derive(Clone, Copy)]
struct ReadMetrics {
    length: usize,
    valid_bases: usize,
    gc_fraction: f64,
    a_fraction: f64,
    c_fraction: f64,
    g_fraction: f64,
    t_fraction: f64,
    base_entropy: f64,
    max_homopolymer_nt: usize,
    mono_dominance: f64,
    di_dominance: f64,
    tri_dominance: f64,
    rc_covered_fraction: f64,
    max_rc_tract_nt: usize,
    rc_tract_count: usize,
    dominant_period_nt: Option<usize>,
    periodicity_score: f64,
    repeat_units: Option<f64>,
}

struct CharacterizationRow {
    fragment_id: u64,
    read_id: String,
    mate: &'static str,
    winners: Vec<u16>,
    selected_winners: Vec<u16>,
    metrics: ReadMetrics,
}

#[derive(Clone, Copy)]
struct SummaryStats {
    n: usize,
    mean: f64,
    sd: f64,
    median: f64,
    q1: f64,
    q3: f64,
    p05: f64,
    p95: f64,
}

const SUMMARY_METRICS: &[&str] = &[
    "length",
    "valid_bases",
    "gc_fraction",
    "a_fraction",
    "c_fraction",
    "g_fraction",
    "t_fraction",
    "base_entropy",
    "max_homopolymer_nt",
    "mono_dominance",
    "di_dominance",
    "tri_dominance",
    "rc_covered_fraction",
    "max_rc_tract_nt",
    "rc_tract_count",
    "dominant_period_nt",
    "periodicity_score",
    "repeat_units",
];

fn characterize_sample(
    assignments_path: PathBuf,
    assignment_map_path: Option<PathBuf>,
    r1_path: PathBuf,
    r2_path: Option<PathBuf>,
    requested_subfamilies: Vec<String>,
    output_prefix: PathBuf,
    include_ambiguous: bool,
    min_rc_tract: usize,
    max_period: usize,
    min_periodicity_score: f64,
    period_equivalence_delta: f64,
    debug: bool,
) {
    assert!(min_rc_tract >= 1, "--min-rc-tract must be at least 1");
    assert!(max_period >= 1, "--max-period must be at least 1");
    assert!(
        (0.0..=1.0).contains(&min_periodicity_score),
        "--min-periodicity-score must be between 0 and 1"
    );
    assert!(
        (0.0..=1.0).contains(&period_equivalence_delta),
        "--period-equivalence-delta must be between 0 and 1"
    );

    let map_path = assignment_map_path
        .unwrap_or_else(|| default_assignment_map_path(&assignments_path));
    let subfamily_names = load_assignment_map(&map_path);

    let mut selected_order: Vec<u16> = Vec::new();
    let mut selected_set: FxHashSet<u16> = FxHashSet::default();
    for requested in requested_subfamilies {
        let match_id = subfamily_names
            .iter()
            .position(|name| name == &requested)
            .unwrap_or_else(|| {
                panic!(
                    "Requested subfamily '{}' is not present in assignment map {:?}",
                    requested, map_path
                )
            });
        assert!(
            match_id <= u16::MAX as usize,
            "Assignment map contains an ID too large for u16"
        );
        let id = match_id as u16;
        if selected_set.insert(id) {
            selected_order.push(id);
        }
    }

    if debug {
        eprintln!(
            "[INFO] Characterization inputs: assignments={:?}, assignment_map={:?}, include_ambiguous={}",
            assignments_path, map_path, include_ambiguous
        );
        eprintln!(
            "[INFO] Characterization settings: min_rc_tract={}, max_period={}, min_periodicity_score={:.3}, period_equivalence_delta={:.3}",
            min_rc_tract,
            max_period,
            min_periodicity_score,
            period_equivalence_delta
        );
    }

    let assignments_file =
        File::open(&assignments_path).expect("Cannot open quant assignment output file");
    let mut assignments_reader =
        BufReader::with_capacity(512 * 1024, assignments_file);
    let mut assignment_line = String::new();

    let per_read_path = append_path_suffix(&output_prefix, ".per_read.tsv");
    let summary_path = append_path_suffix(&output_prefix, ".summary.tsv");
    let mut per_read = BufWriter::new(
        File::create(&per_read_path).expect("Cannot create per-read characterization output"),
    );
    writeln!(
        per_read,
        "fragment_id\tread_id\tmate\tsubfamily\tassignment\tassignment_ambiguous\tlength\tvalid_bases\tgc_fraction\ta_fraction\tc_fraction\tg_fraction\tt_fraction\tbase_entropy\tmax_homopolymer_nt\tmono_dominance\tdi_dominance\ttri_dominance\trc_covered_fraction\tmax_rc_tract_nt\trc_tract_count\tdominant_period_nt\tperiodicity_score\trepeat_units"
    )
    .unwrap();

    let mut summary_rows: FxHashMap<u16, Vec<ReadMetrics>> = FxHashMap::default();
    for &id in &selected_order {
        summary_rows.insert(id, Vec::new());
    }

    let mut r1 = get_fastx_reader(&r1_path);
    let mut r2 = r2_path.as_ref().map(get_fastx_reader);
    let batch_size = 20_000usize;
    let mut fragment_id = 0u64;
    let mut total_characterized_rows = 0u64;

    loop {
        let mut batch = Vec::with_capacity(batch_size);
        let mut records_read = 0usize;

        for _ in 0..batch_size {
            match r1.next() {
                Some(Ok(a)) => {
                    fragment_id += 1;
                    records_read += 1;

                    assignment_line.clear();
                    let assignment_bytes = assignments_reader
                        .read_line(&mut assignment_line)
                        .expect("Failed while reading quant assignment output");
                    if assignment_bytes == 0 {
                        panic!(
                            "Assignment file ended before R1 FASTQ at fragment {}",
                            fragment_id
                        );
                    }
                    let winners = parse_assignment_line(
                        &assignment_line,
                        fragment_id,
                        subfamily_names.len(),
                    );

                    let r2_record = if let Some(r2_reader) = r2.as_mut() {
                        match r2_reader.next() {
                            Some(Ok(b)) => Some((
                                tsv_safe(&String::from_utf8_lossy(b.id())),
                                b.seq().to_vec(),
                            )),
                            Some(Err(_)) => panic!("Failed while reading R2 FASTQ"),
                            None => panic!("R2 FASTQ ended before R1 FASTQ"),
                        }
                    } else {
                        None
                    };

                    // Assignment is a fragment-level product from quant. Characterization
                    // only selects from that product; it never re-runs k-mer assignment.
                    if winners.is_empty() || (!include_ambiguous && winners.len() != 1) {
                        continue;
                    }

                    let selected_winners: Vec<u16> = winners
                        .iter()
                        .copied()
                        .filter(|id| selected_set.contains(id))
                        .collect();
                    if selected_winners.is_empty() {
                        continue;
                    }

                    batch.push(CharacterizationFragment {
                        fragment_id,
                        r1_id: tsv_safe(&String::from_utf8_lossy(a.id())),
                        r1_seq: a.seq().to_vec(),
                        r2: r2_record,
                        winners,
                        selected_winners,
                    });
                }
                Some(Err(_)) => panic!("Failed while reading R1 FASTQ"),
                None => break,
            }
        }

        if records_read == 0 {
            break;
        }

        let grouped_rows: Vec<Vec<CharacterizationRow>> = batch
            .par_iter()
            .map(|fragment| {
                let mut rows = Vec::with_capacity(if fragment.r2.is_some() { 2 } else { 1 });
                rows.push(CharacterizationRow {
                    fragment_id: fragment.fragment_id,
                    read_id: fragment.r1_id.clone(),
                    mate: "R1",
                    winners: fragment.winners.clone(),
                    selected_winners: fragment.selected_winners.clone(),
                    metrics: characterize_read(
                        &fragment.r1_seq,
                        min_rc_tract,
                        max_period,
                        min_periodicity_score,
                        period_equivalence_delta,
                    ),
                });

                if let Some((r2_id, r2_seq)) = &fragment.r2 {
                    rows.push(CharacterizationRow {
                        fragment_id: fragment.fragment_id,
                        read_id: r2_id.clone(),
                        mate: "R2",
                        winners: fragment.winners.clone(),
                        selected_winners: fragment.selected_winners.clone(),
                        metrics: characterize_read(
                            r2_seq,
                            min_rc_tract,
                            max_period,
                            min_periodicity_score,
                            period_equivalence_delta,
                        ),
                    });
                }
                rows
            })
            .collect();

        for rows in grouped_rows {
            for row in rows {
                let selected_names = join_subfamily_names(&subfamily_names, &row.selected_winners);
                let assignment_names = join_subfamily_names(&subfamily_names, &row.winners);
                write_characterization_row(
                    &mut per_read,
                    &row,
                    &selected_names,
                    &assignment_names,
                );
                for selected_id in &row.selected_winners {
                    summary_rows
                        .get_mut(selected_id)
                        .expect("Selected subfamily summary accumulator missing")
                        .push(row.metrics);
                }
                total_characterized_rows += 1;
            }
        }

        if debug && fragment_id % 100_000 < batch_size as u64 {
            eprintln!(
                "[INFO] Processed {} assignment/FASTQ fragments; emitted {} characterized read rows",
                fragment_id, total_characterized_rows
            );
        }
    }

    if let Some(r2_reader) = r2.as_mut() {
        match r2_reader.next() {
            Some(_) => panic!("R2 FASTQ contains more records than R1 FASTQ"),
            None => {}
        }
    }

    assignment_line.clear();
    if assignments_reader
        .read_line(&mut assignment_line)
        .expect("Failed while checking assignment file length")
        != 0
    {
        panic!("Assignment file contains more records than R1 FASTQ");
    }

    per_read.flush().unwrap();
    write_characterization_summary(
        &summary_path,
        &subfamily_names,
        &selected_order,
        &summary_rows,
    );

    println!("[INFO] Per-read characterization: {:?}", per_read_path);
    println!("[INFO] Characterization summary: {:?}", summary_path);
}

fn default_assignment_map_path(assignments_path: &PathBuf) -> PathBuf {
    let mut map_path = assignments_path.clone();
    let file_name = assignments_path
        .file_name()
        .and_then(|name| name.to_str())
        .expect("Assignment path must have a valid UTF-8 file name");
    map_path.set_file_name(format!("{}.map.tsv", file_name));
    map_path
}

fn load_assignment_map(path: &PathBuf) -> Vec<String> {
    let file = File::open(path).unwrap_or_else(|_| {
        panic!(
            "Cannot open assignment map {:?}. Pass --assignment-map explicitly if it is not next to the assignment file.",
            path
        )
    });
    let mut reader = BufReader::new(file);
    let mut line = String::new();

    let n = reader
        .read_line(&mut line)
        .expect("Failed while reading assignment map header");
    assert!(n > 0, "Assignment map is empty");
    assert_eq!(
        line.trim(),
        "ID\tSubfamily",
        "Assignment map must begin with the header 'ID\\tSubfamily'"
    );

    let mut names_by_id: Vec<Option<String>> = Vec::new();
    let mut line_number = 1usize;
    loop {
        line.clear();
        if reader
            .read_line(&mut line)
            .expect("Failed while reading assignment map")
            == 0
        {
            break;
        }
        line_number += 1;
        let trimmed = line.trim_end_matches(&['\r', '\n'][..]);
        if trimmed.is_empty() {
            continue;
        }

        let (id_text, name) = trimmed.split_once('\t').unwrap_or_else(|| {
            panic!(
                "Invalid assignment map row {}: expected ID<TAB>Subfamily",
                line_number
            )
        });
        let id: usize = id_text.parse().unwrap_or_else(|_| {
            panic!(
                "Invalid assignment ID '{}' in map row {}",
                id_text, line_number
            )
        });
        assert!(
            id <= u16::MAX as usize,
            "Assignment ID {} in map row {} exceeds u16 range",
            id,
            line_number
        );
        assert!(
            !name.is_empty(),
            "Empty subfamily name in assignment map row {}",
            line_number
        );

        if names_by_id.len() <= id {
            names_by_id.resize(id + 1, None);
        }
        assert!(
            names_by_id[id].is_none(),
            "Duplicate assignment ID {} in map row {}",
            id,
            line_number
        );
        names_by_id[id] = Some(name.to_string());
    }

    assert!(
        !names_by_id.is_empty(),
        "Assignment map contains no subfamily entries"
    );

    names_by_id
        .into_iter()
        .enumerate()
        .map(|(id, name)| {
            name.unwrap_or_else(|| panic!("Assignment map is missing ID {}", id))
        })
        .collect()
}

fn parse_assignment_line(line: &str, fragment_id: u64, num_subfamilies: usize) -> Vec<u16> {
    let trimmed = line.trim();
    if trimmed == "-1" {
        return Vec::new();
    }
    assert!(
        !trimmed.is_empty(),
        "Empty assignment record at fragment {}",
        fragment_id
    );

    let mut winners = Vec::new();
    for token in trimmed.split(',') {
        let id: u16 = token.parse().unwrap_or_else(|_| {
            panic!(
                "Invalid assignment token '{}' at fragment {}",
                token, fragment_id
            )
        });
        assert!(
            (id as usize) < num_subfamilies,
            "Assignment ID {} at fragment {} is absent from the assignment map",
            id,
            fragment_id
        );
        winners.push(id);
    }
    winners.sort_unstable();
    winners.dedup();
    winners
}

fn append_path_suffix(path: &PathBuf, suffix: &str) -> PathBuf {
    let mut os = path.as_os_str().to_os_string();
    os.push(suffix);
    PathBuf::from(os)
}

fn tsv_safe(text: &str) -> String {
    text.replace('\t', " ").replace('\n', " ").replace('\r', " ")
}

fn join_subfamily_names(names: &[String], ids: &[u16]) -> String {
    ids.iter()
        .map(|&id| names[id as usize].as_str())
        .collect::<Vec<_>>()
        .join(",")
}

fn characterize_read(
    seq: &[u8],
    min_rc_tract: usize,
    max_period: usize,
    min_periodicity_score: f64,
    period_equivalence_delta: f64,
) -> ReadMetrics {
    let length = seq.len();
    let mut counts = [0usize; 4];
    for &base in seq {
        if let Some(code) = base_code(base) {
            counts[code as usize] += 1;
        }
    }
    let valid_bases: usize = counts.iter().sum();
    let denom = valid_bases as f64;
    let a_fraction = safe_fraction(counts[0], valid_bases);
    let c_fraction = safe_fraction(counts[1], valid_bases);
    let g_fraction = safe_fraction(counts[2], valid_bases);
    let t_fraction = safe_fraction(counts[3], valid_bases);
    let gc_fraction = safe_fraction(counts[1] + counts[2], valid_bases);

    let base_entropy = if valid_bases == 0 {
        0.0
    } else {
        let mut entropy = 0.0;
        for &count in &counts {
            if count > 0 {
                let p = count as f64 / denom;
                entropy -= p * p.log2();
            }
        }
        entropy / 2.0
    };

    let max_homopolymer_nt = longest_homopolymer(seq);
    let mono_dominance = dominant_kmer_fraction(seq, 1);
    let di_dominance = dominant_kmer_fraction(seq, 2);
    let tri_dominance = dominant_kmer_fraction(seq, 3);
    let (rc_covered_fraction, max_rc_tract_nt, rc_tract_count) =
        reverse_complementarity_metrics(seq, min_rc_tract, valid_bases);
    let (dominant_period_nt, periodicity_score, repeat_units) = periodicity_metrics(
        seq,
        max_period,
        min_periodicity_score,
        period_equivalence_delta,
    );

    ReadMetrics {
        length,
        valid_bases,
        gc_fraction,
        a_fraction,
        c_fraction,
        g_fraction,
        t_fraction,
        base_entropy,
        max_homopolymer_nt,
        mono_dominance,
        di_dominance,
        tri_dominance,
        rc_covered_fraction,
        max_rc_tract_nt,
        rc_tract_count,
        dominant_period_nt,
        periodicity_score,
        repeat_units,
    }
}

#[inline]
fn base_code(base: u8) -> Option<u8> {
    match base {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

#[inline]
fn safe_fraction(numerator: usize, denominator: usize) -> f64 {
    if denominator == 0 {
        0.0
    } else {
        numerator as f64 / denominator as f64
    }
}

fn longest_homopolymer(seq: &[u8]) -> usize {
    let mut best = 0usize;
    let mut run = 0usize;
    let mut previous: Option<u8> = None;
    for &base in seq {
        match base_code(base) {
            Some(code) if previous == Some(code) => {
                run += 1;
            }
            Some(code) => {
                previous = Some(code);
                run = 1;
            }
            None => {
                previous = None;
                run = 0;
            }
        }
        best = best.max(run);
    }
    best
}

fn dominant_kmer_fraction(seq: &[u8], k: usize) -> f64 {
    if k == 0 || seq.len() < k {
        return 0.0;
    }

    let mut counts: FxHashMap<u8, usize> = FxHashMap::default();
    let mut valid_windows = 0usize;
    for start in 0..=(seq.len() - k) {
        let mut key = 0u8;
        let mut valid = true;
        for &base in &seq[start..start + k] {
            if let Some(code) = base_code(base) {
                key = (key << 2) | code;
            } else {
                valid = false;
                break;
            }
        }
        if valid {
            *counts.entry(key).or_insert(0) += 1;
            valid_windows += 1;
        }
    }

    if valid_windows == 0 {
        0.0
    } else {
        let max_count = counts.values().copied().max().unwrap_or(0);
        max_count as f64 / valid_windows as f64
    }
}

fn reverse_complementarity_metrics(
    seq: &[u8],
    min_tract: usize,
    valid_bases: usize,
) -> (f64, usize, usize) {
    let n = seq.len();
    if n < 2 || min_tract == 0 {
        return (0.0, 0, 0);
    }

    let mut covered = vec![false; n];
    let mut max_tract = 0usize;
    let mut tract_count = 0usize;

    // A reverse-complement tract is a maximal run along an anti-diagonal:
    // (i, j), (i+1, j-1), ... where the bases are exact complements.
    // Scanning maximal runs prevents a long tract from being counted as many
    // overlapping short k-mer matches.
    for sum in 1..=(2 * n - 3) {
        let first_i = sum.saturating_sub(n - 1);
        let last_i = std::cmp::min((sum - 1) / 2, n - 1);
        if first_i > last_i {
            continue;
        }

        let mut run_start = 0usize;
        let mut run_len = 0usize;
        for i in first_i..=last_i {
            let j = sum - i;
            let complementary = match (base_code(seq[i]), base_code(seq[j])) {
                (Some(a), Some(b)) => (a ^ b) == 3,
                _ => false,
            };

            if complementary {
                if run_len == 0 {
                    run_start = i;
                }
                run_len += 1;
            } else if run_len > 0 {
                finalize_rc_run(
                    sum,
                    run_start,
                    run_len,
                    min_tract,
                    &mut covered,
                    &mut max_tract,
                    &mut tract_count,
                );
                run_len = 0;
            }
        }

        if run_len > 0 {
            finalize_rc_run(
                sum,
                run_start,
                run_len,
                min_tract,
                &mut covered,
                &mut max_tract,
                &mut tract_count,
            );
        }
    }

    let covered_bases = covered.iter().filter(|&&is_covered| is_covered).count();
    (
        safe_fraction(covered_bases, valid_bases),
        max_tract,
        tract_count,
    )
}

fn finalize_rc_run(
    sum: usize,
    run_start: usize,
    run_len: usize,
    min_tract: usize,
    covered: &mut [bool],
    max_tract: &mut usize,
    tract_count: &mut usize,
) {
    if run_len < min_tract {
        return;
    }
    *tract_count += 1;
    *max_tract = (*max_tract).max(run_len);
    for offset in 0..run_len {
        let i = run_start + offset;
        let j = sum - i;
        covered[i] = true;
        covered[j] = true;
    }
}

fn periodicity_metrics(
    seq: &[u8],
    max_period: usize,
    min_score: f64,
    equivalence_delta: f64,
) -> (Option<usize>, f64, Option<f64>) {
    let n = seq.len();
    let max_test_period = std::cmp::min(max_period, n / 2);
    if max_test_period == 0 {
        return (None, 0.0, None);
    }

    let mut scores = Vec::with_capacity(max_test_period);
    let mut best_score = 0.0f64;
    for period in 1..=max_test_period {
        let mut comparable = 0usize;
        let mut matches = 0usize;
        for i in period..n {
            if let (Some(a), Some(b)) = (base_code(seq[i]), base_code(seq[i - period])) {
                comparable += 1;
                if a == b {
                    matches += 1;
                }
            }
        }
        let score = safe_fraction(matches, comparable);
        best_score = best_score.max(score);
        scores.push(score);
    }

    if best_score < min_score {
        return (None, best_score, None);
    }

    let cutoff = (best_score - equivalence_delta).max(min_score);
    let dominant_period = scores
        .iter()
        .enumerate()
        .find_map(|(idx, &score)| if score >= cutoff { Some(idx + 1) } else { None });

    match dominant_period {
        Some(period) => (
            Some(period),
            scores[period - 1],
            Some(n as f64 / period as f64),
        ),
        None => (None, best_score, None),
    }
}

fn write_characterization_row(
    writer: &mut BufWriter<File>,
    row: &CharacterizationRow,
    selected_names: &str,
    assignment_names: &str,
) {
    let m = row.metrics;
    let period = m
        .dominant_period_nt
        .map(|p| p.to_string())
        .unwrap_or_else(|| "NA".to_string());
    let repeat_units = m
        .repeat_units
        .map(|v| format!("{:.6}", v))
        .unwrap_or_else(|| "NA".to_string());

    writeln!(
        writer,
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{}\t{}\t{}\t{:.6}\t{}",
        row.fragment_id,
        row.read_id,
        row.mate,
        selected_names,
        assignment_names,
        row.winners.len() > 1,
        m.length,
        m.valid_bases,
        m.gc_fraction,
        m.a_fraction,
        m.c_fraction,
        m.g_fraction,
        m.t_fraction,
        m.base_entropy,
        m.max_homopolymer_nt,
        m.mono_dominance,
        m.di_dominance,
        m.tri_dominance,
        m.rc_covered_fraction,
        m.max_rc_tract_nt,
        m.rc_tract_count,
        period,
        m.periodicity_score,
        repeat_units
    )
    .unwrap();
}

fn metric_value(metrics: &ReadMetrics, metric: &str) -> Option<f64> {
    match metric {
        "length" => Some(metrics.length as f64),
        "valid_bases" => Some(metrics.valid_bases as f64),
        "gc_fraction" => Some(metrics.gc_fraction),
        "a_fraction" => Some(metrics.a_fraction),
        "c_fraction" => Some(metrics.c_fraction),
        "g_fraction" => Some(metrics.g_fraction),
        "t_fraction" => Some(metrics.t_fraction),
        "base_entropy" => Some(metrics.base_entropy),
        "max_homopolymer_nt" => Some(metrics.max_homopolymer_nt as f64),
        "mono_dominance" => Some(metrics.mono_dominance),
        "di_dominance" => Some(metrics.di_dominance),
        "tri_dominance" => Some(metrics.tri_dominance),
        "rc_covered_fraction" => Some(metrics.rc_covered_fraction),
        "max_rc_tract_nt" => Some(metrics.max_rc_tract_nt as f64),
        "rc_tract_count" => Some(metrics.rc_tract_count as f64),
        "dominant_period_nt" => metrics.dominant_period_nt.map(|v| v as f64),
        "periodicity_score" => Some(metrics.periodicity_score),
        "repeat_units" => metrics.repeat_units,
        _ => None,
    }
}

fn write_characterization_summary(
    path: &PathBuf,
    names: &[String],
    selected_order: &[u16],
    rows_by_subfamily: &FxHashMap<u16, Vec<ReadMetrics>>,
) {
    let mut out = BufWriter::new(File::create(path).expect("Cannot create summary output"));
    writeln!(
        out,
        "subfamily\tmetric\tn\tmean\tsd\tmedian\tq1\tq3\tp05\tp95"
    )
    .unwrap();

    for &id in selected_order {
        let name = names[id as usize].as_str();
        let rows = rows_by_subfamily
            .get(&id)
            .expect("Selected subfamily summary rows missing");
        for &metric in SUMMARY_METRICS {
            let values: Vec<f64> = rows
                .iter()
                .filter_map(|m| metric_value(m, metric))
                .collect();
            if let Some(stats) = summarize_values(values) {
                writeln!(
                    out,
                    "{}\t{}\t{}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}",
                    name,
                    metric,
                    stats.n,
                    stats.mean,
                    stats.sd,
                    stats.median,
                    stats.q1,
                    stats.q3,
                    stats.p05,
                    stats.p95
                )
                .unwrap();
            } else {
                writeln!(
                    out,
                    "{}\t{}\t0\tNA\tNA\tNA\tNA\tNA\tNA\tNA",
                    name, metric
                )
                .unwrap();
            }
        }
    }
    out.flush().unwrap();
}

fn summarize_values(mut values: Vec<f64>) -> Option<SummaryStats> {
    if values.is_empty() {
        return None;
    }
    values.sort_unstable_by(|a, b| a.total_cmp(b));
    let n = values.len();
    let mean = values.iter().sum::<f64>() / n as f64;
    let variance = values
        .iter()
        .map(|v| {
            let d = v - mean;
            d * d
        })
        .sum::<f64>()
        / n as f64;

    Some(SummaryStats {
        n,
        mean,
        sd: variance.sqrt(),
        median: quantile_sorted(&values, 0.50),
        q1: quantile_sorted(&values, 0.25),
        q3: quantile_sorted(&values, 0.75),
        p05: quantile_sorted(&values, 0.05),
        p95: quantile_sorted(&values, 0.95),
    })
}

fn quantile_sorted(values: &[f64], q: f64) -> f64 {
    if values.len() == 1 {
        return values[0];
    }
    let position = q * (values.len() - 1) as f64;
    let lower = position.floor() as usize;
    let upper = position.ceil() as usize;
    if lower == upper {
        values[lower]
    } else {
        let weight = position - lower as f64;
        values[lower] * (1.0 - weight) + values[upper] * weight
    }
}

fn process_batch_rkyv(
    batch: &[(Vec<u8>, Option<Vec<u8>>)],
    k: usize,
    partitions: &[&rkyv::collections::ArchivedHashMap<u128, rkyv::vec::ArchivedVec<u16>>],
    num_classes: usize,
    min_matches: u32,
    track_assignments: bool,
) -> (Vec<f64>, Option<Vec<String>>) {
    let results: Vec<(Vec<(usize, f64)>, String)> = batch
        .par_iter()
        .map(|(s1, s2)| {
            let mut votes = FxHashMap::default();
            scan_read(s1, k, partitions, &mut votes);
            if let Some(s) = s2 {
                scan_read(s, k, partitions, &mut votes);
            }

            if votes.is_empty() {
                (Vec::new(), "-1".to_string())
            } else {
                let max = *votes.values().max().unwrap();
                if max >= min_matches {
                    let winners: Vec<_> = votes
                        .iter()
                        .filter(|&(_, &c)| c == max)
                        .map(|(&i, _)| i)
                        .collect();
                    let inc = 1.0 / winners.len() as f64;
                    let increments: Vec<(usize, f64)> =
                        winners.iter().map(|&id| (id as usize, inc)).collect();

                    let assign_str = if track_assignments {
                        winners
                            .iter()
                            .map(|id| id.to_string())
                            .collect::<Vec<_>>()
                            .join(",")
                    } else {
                        String::new()
                    };

                    (increments, assign_str)
                } else {
                    (Vec::new(), "-1".to_string())
                }
            }
        })
        .collect();

    let mut batch_counts = vec![0.0; num_classes];
    let mut batch_assignments = if track_assignments {
        Some(Vec::with_capacity(batch.len()))
    } else {
        None
    };

    for (increments, assign_str) in results {
        for (id, inc) in increments {
            batch_counts[id] += inc;
        }
        if let Some(assignments) = &mut batch_assignments {
            assignments.push(assign_str);
        }
    }

    (batch_counts, batch_assignments)
}

fn scan_read(
    seq: &[u8],
    k: usize,
    partitions: &[&rkyv::collections::ArchivedHashMap<u128, rkyv::vec::ArchivedVec<u16>>],
    votes: &mut FxHashMap<u16, u32>,
) {
    if seq.len() < k {
        return;
    }

    for i in 0..=(seq.len() - k) {
        if let Some(forward_enc) = encode_kmer(&seq[i..i + k]) {
            let rc_enc = get_rc_kmer(forward_enc, k);
            let canonical = std::cmp::min(forward_enc, rc_enc);
            let part_idx = (canonical & 0xFF) as usize;

            if let Some(ids) = partitions[part_idx].get(&canonical) {
                for &id in ids.as_slice() {
                    *votes.entry(id).or_insert(0) += 1;
                }
            }
        }
    }
}

fn encode_kmer(seq: &[u8]) -> Option<u128> {
    let mut val: u128 = 0;
    for &base in seq {
        val <<= 2;
        match base {
            b'A' | b'a' => val |= 0,
            b'C' | b'c' => val |= 1,
            b'G' | b'g' => val |= 2,
            b'T' | b't' => val |= 3,
            _ => return None,
        }
    }
    Some(val)
}

#[inline(always)]
fn get_rc_kmer(val: u128, k: usize) -> u128 {
    let rev = val.reverse_bits();
    let swapped = ((rev & 0xAAAA_AAAA_AAAA_AAAA_AAAA_AAAA_AAAA_AAAA) >> 1)
        | ((rev & 0x5555_5555_5555_5555_5555_5555_5555_5555) << 1);
    (!swapped) >> (128 - 2 * k)
}


#[cfg(test)]
mod characterization_tests {
    use super::*;

    fn approx_eq(a: f64, b: f64) {
        assert!((a - b).abs() < 1e-9, "{} != {}", a, b);
    }

    #[test]
    fn composition_ignores_ambiguous_bases() {
        let m = characterize_read(b"ACGTN", 5, 12, 0.80, 0.02);
        assert_eq!(m.length, 5);
        assert_eq!(m.valid_bases, 4);
        approx_eq(m.gc_fraction, 0.5);
        approx_eq(m.a_fraction, 0.25);
        approx_eq(m.c_fraction, 0.25);
        approx_eq(m.g_fraction, 0.25);
        approx_eq(m.t_fraction, 0.25);
        approx_eq(m.base_entropy, 1.0);
    }

    #[test]
    fn ambiguous_base_breaks_homopolymer() {
        assert_eq!(longest_homopolymer(b"AAAANAAAA"), 4);
    }

    #[test]
    fn reverse_complement_tract_is_merged_and_covered() {
        let (fraction, max_tract, count) =
            reverse_complementarity_metrics(b"ACGTTGNNNNCAACGT", 6, 12);
        approx_eq(fraction, 1.0);
        assert_eq!(max_tract, 6);
        assert_eq!(count, 1);
    }

    #[test]
    fn primitive_period_is_preferred() {
        let (period, score, units) = periodicity_metrics(b"ATATATAT", 12, 0.80, 0.02);
        assert_eq!(period, Some(2));
        approx_eq(score, 1.0);
        approx_eq(units.unwrap(), 4.0);
    }

    #[test]
    fn example_periods_are_detected() {
        assert_eq!(periodicity_metrics(b"AAAAAAAA", 12, 0.80, 0.02).0, Some(1));
        assert_eq!(periodicity_metrics(b"CAGCAGCAG", 12, 0.80, 0.02).0, Some(3));
        assert_eq!(
            periodicity_metrics(b"AATAACAATAAC", 12, 0.80, 0.02).0,
            Some(6)
        );
    }

    #[test]
    fn exported_assignment_lines_are_parsed() {
        assert!(parse_assignment_line("-1\n", 1, 10).is_empty());
        assert_eq!(parse_assignment_line("3\n", 2, 10), vec![3]);
        assert_eq!(parse_assignment_line("7,2\n", 3, 10), vec![2, 7]);
    }

    #[test]
    fn assignment_map_path_matches_quant_output_name() {
        let assignments = PathBuf::from("results/sample.assignments.tsv");
        assert_eq!(
            default_assignment_map_path(&assignments),
            PathBuf::from("results/sample.assignments.tsv.map.tsv")
        );
    }
}

fn get_fastx_reader(path: &PathBuf) -> Box<dyn FastxReader> {
    let is_gzipped = path
        .extension()
        .and_then(|s| s.to_str())
        .map(|s| s.eq_ignore_ascii_case("gz"))
        .unwrap_or(false);

    if is_gzipped {
        let child = Command::new("unpigz")
            .arg("-c")
            .arg(path)
            .stdout(Stdio::piped())
            .spawn()
            .expect("unpigz missing");
        parse_fastx_reader(BufReader::with_capacity(512 * 1024, child.stdout.unwrap())).unwrap()
    } else {
        needletail::parse_fastx_file(path).unwrap()
    }
}

