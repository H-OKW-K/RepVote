use clap::{Parser, Subcommand};
use needletail::{parse_fastx_reader, FastxReader};
use rayon::prelude::*;
use rkyv::{Archive, Deserialize, Serialize};
use memmap2::MmapOptions;
use rustc_hash::FxHashMap;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Write};
use std::path::PathBuf;
use std::process::{Command, Stdio};
use std::sync::mpsc::sync_channel;
use std::thread;

#[derive(Parser)]
#[command(name = "repeat_counter")]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    Index {
        #[arg(long)]
        file_list: PathBuf,
        #[arg(long)]
        output: PathBuf,
        #[arg(short, long)]
        k: usize,
    },
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
        #[arg(long)]
        debug: bool,
    },
}

// Separate data structures for Partitioning
#[derive(Archive, Deserialize, Serialize)]
#[archive(check_bytes)]
pub struct Partition {
    pub map: HashMap<u64, Vec<u16>>,
}

#[derive(Archive, Deserialize, Serialize)]
#[archive(check_bytes)]
pub struct Metadata {
    pub k: u8,
    pub names: Vec<String>,
    pub offsets: Vec<(u64, u64)>, // Start and end bytes of each of the 256 partitions
}

fn main() {
    let args = Cli::parse();

    match args.command {
        Commands::Index { file_list, output, k } => {
            let k_u8 = u8::try_from(k).expect("k must fit in u8");
            build_index(file_list, output, k_u8);
        }
        Commands::Quant { index, r1, r2, output, threads, debug } => {
            rayon::ThreadPoolBuilder::new().num_threads(threads).build_global().unwrap();
            quantify_sample(index, r1, r2, output, debug);
        }
    }
}

fn build_index(list_path: PathBuf, output: PathBuf, k: u8) {
    println!("[INFO] Starting Disk-Spilling Extraction Phase...");
    let mut names: Vec<String> = Vec::new();

    // Create a temporary directory for the 256 partitions
    let mut temp_dir = output.clone();
    temp_dir.set_extension("tmp_partitions");
    std::fs::create_dir_all(&temp_dir).expect("Failed to create temp directory");

    // Open 256 BufWriters pointing to temp files
    let mut writers: Vec<std::io::BufWriter<File>> = (0..256)
        .map(|i| {
            let p = temp_dir.join(format!("part_{:03}.bin", i));
            std::io::BufWriter::new(File::create(p).unwrap())
        })
        .collect();

    let reader = BufReader::new(File::open(&list_path).unwrap());
    let mut valid_id_counter: u16 = 0;

    // 1. Extract & Route to Disk
    for l in reader.lines() {
        let path = PathBuf::from(l.unwrap());
        
        let mut r = match needletail::parse_fastx_file(&path) {
            Ok(r) => r,
            Err(_) => {
                eprintln!("Skipping invalid/empty file: {:?}", path);
                continue;
            }
        };

        names.push(path.file_stem().unwrap().to_str().unwrap().to_string());
        let current_id = valid_id_counter;
        valid_id_counter += 1;

        while let Some(Ok(rec)) = r.next() {
            let seq = rec.seq();
            if seq.len() == k as usize {
                // Input is already canonicalized unique k-mers, no reverse-complement needed here
                if let Some(enc) = encode_kmer(&seq) {
                    // Extract lowest 8 bits to route perfectly to 256 buckets regardless of K
                    let part_idx = (enc & 0xFF) as usize;
                    
                    // Write raw binary to temp file
                    writers[part_idx].write_all(&enc.to_le_bytes()).unwrap();
                    writers[part_idx].write_all(&current_id.to_le_bytes()).unwrap();
                }
            }
        }
    }

    // Flush and close all temporary files safely
    drop(writers);

    println!("[INFO] Extraction complete. Merging partitions into final index...");

    // 2. Merge & Serialize Phase
    let mut out_writer = std::io::BufWriter::new(File::create(&output).unwrap());
    let mut offsets = Vec::with_capacity(256);
    let mut current_offset: u64 = 0;

    for i in 0..256 {
        let p = temp_dir.join(format!("part_{:03}.bin", i));
        let mut f = File::open(&p).unwrap();
        let mut map: HashMap<u64, Vec<u16>> = HashMap::new();

        // Read the binary stream in 10-byte chunks (8 bytes u64, 2 bytes u16)
        let mut buf = [0u8; 10];
        while f.read_exact(&mut buf).is_ok() {
            let kmer = u64::from_le_bytes(buf[0..8].try_into().unwrap());
            let id = u16::from_le_bytes(buf[8..10].try_into().unwrap());
            map.entry(kmer).or_default().push(id);
        }

        // Serialize just this partition chunk
        let partition = Partition { map };
        let bytes = rkyv::to_bytes::<_, 4096>(&partition).unwrap();
        
        // Write out and track offsets
        out_writer.write_all(&bytes).unwrap();
        let chunk_len = bytes.len() as u64;
        offsets.push((current_offset, current_offset + chunk_len));
        current_offset += chunk_len;

        // Cleanup the processed temp file to keep disk usage low
        std::fs::remove_file(p).unwrap();
    }

    // 3. Write Metadata and Footer Offsets
    let meta = Metadata { k, names, offsets };
    let meta_bytes = rkyv::to_bytes::<_, 4096>(&meta).unwrap();
    out_writer.write_all(&meta_bytes).unwrap();
    
    // Write 8 byte pointer to tell the reader exactly where the metadata archive starts
    let meta_offset = current_offset; 
    out_writer.write_all(&meta_offset.to_le_bytes()).unwrap();

    // Clean up empty directory
    std::fs::remove_dir(temp_dir).unwrap();
    println!("[INFO] Index successfully written to Zero-Copy format.");
}

fn quantify_sample(
    index_path: PathBuf,
    r1_path: PathBuf,
    r2_path: Option<PathBuf>,
    output_path: PathBuf,
    debug: bool,
) {
    if debug { eprintln!("[INFO] Mapping index..."); }

    // 1. Instant Load (Zero-Copy)
    let file = File::open(&index_path).expect("Index not found");
    let mmap = unsafe { MmapOptions::new().map(&file).expect("Failed to mmap") };
    let total_len = mmap.len();

    // Read the last 8 bytes to find where our metadata chunk starts
    let meta_ptr_bytes: [u8; 8] = mmap[total_len - 8..total_len].try_into().unwrap();
    let meta_offset = u64::from_le_bytes(meta_ptr_bytes) as usize;

    // Load Metadata
    let meta = unsafe { rkyv::archived_root::<Metadata>(&mmap[meta_offset..total_len - 8]) };
    let k = meta.k as usize;
    let num_classes = meta.names.len();

    // Construct the Partition pointers array
    let mut partitions = Vec::with_capacity(256);
    for i in 0..256 {
        let start = meta.offsets[i].0 as usize;
        let end = meta.offsets[i].1 as usize;
        let archived_part = unsafe { rkyv::archived_root::<Partition>(&mmap[start..end]) };
        partitions.push(&archived_part.map);
    }

    let mut final_counts = vec![0.0; num_classes];

    // Threaded Reader -> Channel -> Rayon Compute
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
                if batch.is_empty() { break; }
                if tx.send(batch).is_err() { break; }
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
                if batch.is_empty() { break; }
                if tx.send(batch).is_err() { break; }
            }
        }
    });

    for batch in rx {
        let batch_counts = process_batch_rkyv(&batch, k, &partitions, num_classes);
        for (i, c) in batch_counts.iter().enumerate() {
            final_counts[i] += c;
        }
    }

    if debug { eprintln!("[INFO] Finished processing."); }

    let mut out = File::create(output_path).unwrap();
    writeln!(out, "RepeatType\tCount").unwrap();
    for (i, count) in final_counts.iter().enumerate() {
        writeln!(out, "{}\t{:.2}", meta.names[i], count).unwrap();
    }
}

fn process_batch_rkyv(
    batch: &[(Vec<u8>, Option<Vec<u8>>)],
    k: usize,
    partitions: &[&rkyv::collections::ArchivedHashMap<u64, rkyv::vec::ArchivedVec<u16>>],
    num_classes: usize,
) -> Vec<f64> {
    batch
        .par_iter()
        .fold(
            || vec![0.0; num_classes],
            |mut acc, (s1, s2)| {
                let mut votes = FxHashMap::default();
                scan_read(s1, k, partitions, &mut votes);
                if let Some(s) = s2 {
                    scan_read(s, k, partitions, &mut votes);
                }

                if !votes.is_empty() {
                    let max = *votes.values().max().unwrap();
                    let winners: Vec<_> = votes
                        .iter()
                        .filter(|&(_, &c)| c == max)
                        .map(|(&i, _)| i)
                        .collect();
                    let inc = 1.0 / winners.len() as f64;
                    for id in winners {
                        acc[id as usize] += inc;
                    }
                }
                acc
            },
        )
        .reduce(
            || vec![0.0; num_classes],
            |mut a, b| {
                for (i, x) in b.iter().enumerate() {
                    a[i] += x;
                }
                a
            },
        )
}

fn scan_read(
    seq: &[u8],
    k: usize,
    partitions: &[&rkyv::collections::ArchivedHashMap<u64, rkyv::vec::ArchivedVec<u16>>],
    votes: &mut FxHashMap<u16, u32>,
) {
    if seq.len() < k { return; }

    for i in 0..=(seq.len() - k) {
        if let Some(forward_enc) = encode_kmer(&seq[i..i + k]) {
            let rc_enc = get_rc_kmer(forward_enc, k);
            let canonical = std::cmp::min(forward_enc, rc_enc);

            // Find correct partition slice using the lowest 8 bits
            let part_idx = (canonical & 0xFF) as usize;

            if let Some(ids) = partitions[part_idx].get(&canonical) {
                for &id in ids.as_slice() {
                    *votes.entry(id).or_insert(0) += 1;
                }
            }
        }
    }
}

fn encode_kmer(seq: &[u8]) -> Option<u64> {
    let mut val: u64 = 0;
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

fn get_rc_kmer(mut val: u64, k: usize) -> u64 {
    let mut rc = 0;
    for _ in 0..k {
        let base = val & 3;
        rc = (rc << 2) | ((!base) & 3);
        val >>= 2;
    }
    rc
}
