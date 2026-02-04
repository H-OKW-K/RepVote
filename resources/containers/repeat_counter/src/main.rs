use clap::{Parser, Subcommand};
use needletail::{parse_fastx_reader, FastxReader};
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::PathBuf;
use rayon::prelude::*;
use std::time::Instant;
use std::sync::mpsc::sync_channel;
use std::thread;
use std::process::{Command, Stdio};
use std::collections::HashMap;

// rkyv imports
use rkyv::{Archive, Deserialize, Serialize};
use memmap2::MmapOptions;

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
    // RESTORED: "Quant" mode for single-sample processing (matches Snakemake)
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

#[derive(Archive, Deserialize, Serialize)]
#[archive(check_bytes)]
struct Database {
    k: u8,
    names: Vec<String>,
    map: HashMap<u64, Vec<u16>>, 
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

            // 1. Instant Load (Zero-Copy)
            if debug { eprintln!("[INFO] Mapping index..."); }
            let file = File::open(&index).expect("Index not found");
            let mmap = unsafe { MmapOptions::new().map(&file).expect("Failed to mmap") };
            let db = unsafe { rkyv::archived_root::<Database>(&mmap) };

            // 2. Process Single Sample
            quantify_sample(db, r1, r2, output, debug);
        }
    }
}

fn build_index(list_path: PathBuf, output: PathBuf, k: u8) {
    println!("Building index...");
    let mut map: HashMap<u64, Vec<u16>> = HashMap::new();
    let mut names: Vec<String> = Vec::new();

    let reader = BufReader::new(File::open(&list_path).unwrap());
    
    // --- FIX: Use a manual counter, not the line number ---
    let mut valid_id_counter: u16 = 0; 

    for l in reader.lines() {
        let path = PathBuf::from(l.unwrap());
        
        // 1. Try to open the file
        let mut r = match needletail::parse_fastx_file(&path) {
            Ok(r) => r,
            Err(_) => {
                eprintln!("Skipping invalid/empty file: {:?}", path);
                continue; // Do NOT increment the counter
            }, 
        };
        
        // 2. Only now do we register the name and ID
        names.push(path.file_stem().unwrap().to_str().unwrap().to_string());
        let current_id = valid_id_counter;
        valid_id_counter += 1; // Increment for the next valid file

        while let Some(Ok(rec)) = r.next() {
            let seq = rec.seq();
            if seq.len() == k as usize {
                if let Some(enc) = encode_kmer(&seq) {
                    map.entry(enc).or_default().push(current_id);
                }
            }
        }
    }

    println!("Saving Zero-Copy Index...");
    let mut writer = File::create(output).unwrap();
    let bytes = rkyv::to_bytes::<_, 4096>(&Database { k, names, map }).unwrap();
    writer.write_all(&bytes).unwrap();
    println!("Done.");
}

fn quantify_sample(
    db: &rkyv::Archived<Database>, 
    r1_path: PathBuf, 
    r2_path: Option<PathBuf>, 
    output_path: PathBuf,
    debug: bool
) {
    let k = db.k as usize;
    let num_classes = db.names.len();
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
                        (Some(Ok(a)), Some(Ok(b))) => batch.push((a.seq().to_vec(), Some(b.seq().to_vec()))),
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
                    if let Some(Ok(a)) = r1.next() { batch.push((a.seq().to_vec(), None)); } else { break; }
                }
                if batch.is_empty() { break; }
                if tx.send(batch).is_err() { break; }
            }
        }
    });

    for batch in rx {
        let batch_counts = process_batch_rkyv(&batch, k, &db.map, num_classes);
        for (i, c) in batch_counts.iter().enumerate() {
            final_counts[i] += c;
        }
    }

    if debug { eprintln!("[INFO] Finished processing."); }

    let mut out = File::create(output_path).unwrap();
    writeln!(out, "RepeatType\tCount").unwrap();
    for (i, count) in final_counts.iter().enumerate() {
        writeln!(out, "{}\t{:.2}", db.names[i], count).unwrap();
    }
}

fn process_batch_rkyv(
    batch: &Vec<(Vec<u8>, Option<Vec<u8>>)>, 
    k: usize, 
    map: &rkyv::collections::ArchivedHashMap<u64, rkyv::vec::ArchivedVec<u16>>, 
    num_classes: usize
) -> Vec<f64> {
    batch.par_iter().fold(
        || vec![0.0; num_classes],
        |mut acc, (s1, s2)| {
            let mut votes = rustc_hash::FxHashMap::default();
            scan_read(s1, k, map, &mut votes);
            if let Some(s) = s2 { scan_read(s, k, map, &mut votes); }

            if !votes.is_empty() {
                let max = *votes.values().max().unwrap();
                let winners: Vec<_> = votes.iter().filter(|&(_, &c)| c == max).map(|(&i, _)| i).collect();
                let inc = 1.0 / winners.len() as f64;
                for id in winners { acc[id as usize] += inc; }
            }
            acc
        }
    ).reduce(
        || vec![0.0; num_classes], 
        |mut a, b| { for (i, x) in b.iter().enumerate() { a[i] += x; } a }
    )
}

fn scan_read(
    seq: &[u8], 
    k: usize, 
    map: &rkyv::collections::ArchivedHashMap<u64, rkyv::vec::ArchivedVec<u16>>, 
    votes: &mut rustc_hash::FxHashMap<u16, u32>
) {
    if seq.len() < k { return; }
    
    // Iterate through every valid window in the read
    for i in 0..=(seq.len() - k) {
        if let Some(forward_enc) = encode_kmer(&seq[i..i+k]) {
            // 1. Calculate Reverse Complement integer
            let rc_enc = get_rc_kmer(forward_enc, k);
            
            // 2. Canonicalize: choose the smaller integer (lexicographically first)
            // This ensures matches regardless of which strand was sequenced
            let canonical = std::cmp::min(forward_enc, rc_enc);

            // 3. Lookup
            if let Some(ids) = map.get(&canonical) {
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
            b'A'|b'a'=> val |=0, b'C'|b'c'=> val |=1, b'G'|b'g'=> val |=2, b'T'|b't'=> val |=3, _=> return None,
        }
    }
    Some(val)
}

fn get_fastx_reader(path: &PathBuf) -> Box<dyn FastxReader> {
    let is_gzipped = path.extension().and_then(|s| s.to_str()).map(|s| s.eq_ignore_ascii_case("gz")).unwrap_or(false);
    if is_gzipped {
        let child = Command::new("unpigz").arg("-c").arg(path).stdout(Stdio::piped()).spawn().expect("unpigz missing");
        parse_fastx_reader(BufReader::with_capacity(512*1024, child.stdout.unwrap())).unwrap()
    } else {
        needletail::parse_fastx_file(path).unwrap()
    }
}

fn get_rc_kmer(mut val: u64, k: usize) -> u64 {
    let mut rc = 0;
    for _ in 0..k {
        let base = val & 3;      // Get last base (lowest 2 bits)
        rc = (rc << 2) | ((!base) & 3); // Complement it and push to accumulator
        val >>= 2;               // Drop processed base
    }
    rc
}