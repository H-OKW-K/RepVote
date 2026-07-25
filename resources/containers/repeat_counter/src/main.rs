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
        repeat_intervals_by_chrom
            .entry(chrom.clone())
            .or_default()
            .push((start, end));

        let family = get_family_from_name(raw_name);
        if family == "Unknown" || family == "Simple_repeat" || family == "Low_complexity" {
            continue;
        }

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
        let mut map_path = path.clone();
        let file_name = format!(
            "{}.map.tsv",
            path.file_name().unwrap().to_str().unwrap()
        );
        map_path.set_file_name(file_name);

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
