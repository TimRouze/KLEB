#![allow(dead_code)]

mod bloom;
mod utils;
mod kmer;
mod lock;
mod parallel_bloom;
use clap::Parser;
use rand::random;
use seq_io::fasta::Reader;
use seq_io::BaseRecord;
use std::collections::HashSet;
use std::sync::{Arc, Mutex, RwLock};
use std::error::Error;
use std::fs::File;
use std::io::{self, BufRead};
use std::cmp::min;//bebou
use std::env;
use::csv::Writer;
use::rayon::prelude::*;
use bloom::{BloomFilter, AggregatingBloomFilter};
use kmer::{Kmer, RawKmer, Base};

use crate::utils::num2str;
//use parallel_bloom::AggregatingBloomFilter;

const K: usize = 31;
const BLOCK_SIZE: usize = 1 << (12 - 3);
const SHARD_AMOUNT: usize = 1024;
const KMER_BITS: usize = 2*K;
pub type KT = u64;

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Input file of file (.txt)
    input: String,
    /// Output file (defaults to out.csv)
    #[arg(short, long, default_value_t = String::from("out.csv"))]
    output: String,
    /// Number of threads (defaults to all available threads)
    #[arg(short, long, default_value_t = 1)]
    threads: usize,
    /// Memory (in GB) allocated to Bloom filters (defaults to 4GB)
    #[arg(short, long, default_value_t = 4)]
    memory: usize,
    /// Number of hashes used in Bloom filters
    #[arg(short = 'H', long, default_value_t = 1)]
    hashes: usize,
    /// Seed used for hash functions
    #[arg(short, long, default_value_t = 101010)]
    seed: u64,
    /// Modimizer ?
    #[arg(short = 'M', long, default_value_t = 1)]
    modimizer: u64,
    /// Number of Agregated Bloom Filter
    #[arg(short, long, default_value_t = 2)]
    filters: u64,
}
//TODO OPTION NB ABF.
//TODO add K as args ?
fn main() {
    let args = Args::parse();
    let input_fof = args.input.as_str();
    let size =  args.memory * 1_000_000_000;
    let hashes = args.hashes;
    let seed = args.seed as u64;
    let nb_filters = args.filters as u64;
    env::set_var("RAYON_NUM_THREADS", args.threads.to_string());
    if let Ok(lines_temp) = utils::read_lines(input_fof){
        let nb_files = lines_temp.count();
        match process_fof_parallel(input_fof, args.modimizer, nb_files, size, hashes, seed, nb_filters) {
            Ok(hist_mutex) => {
                println!("All {} files have been read...\nWriting output...", nb_files);
                write_output(hist_mutex, nb_files, args.output).unwrap();
            }
            Err(err) => eprintln!("Error reading or processing file: {}", err),
        }
    }
}
fn process_fof_parallel(filename: &str, modimizer: u64, nb_files: usize, size: usize, hashes: usize, seed: u64, nb_filters: u64) -> io::Result<Vec<Arc<Mutex<u64>>>>{
    let file = File::open(filename)?;
    let reader = io::BufReader::new(file);
    //let mut mutex_vect = Vec::new();
    /*for i in 1..nb_filters{
        mutex_vect.push(Arc::new(Mutex::new(AggregatingBloomFilter::new_with_seed(size, hashes, seed+random::<u64>()))));
    }*/
    let mut agregated_bf_vector: Vec<Arc<Mutex<AggregatingBloomFilter>>> = Vec::new();
    let mut agregated_bf_vector_2: Vec<Arc<Mutex<AggregatingBloomFilter>>> = Vec::new();
    let shard_size = size/SHARD_AMOUNT;
    for _ in 0..SHARD_AMOUNT{
        agregated_bf_vector.push(Arc::new(Mutex::new(AggregatingBloomFilter::new_with_seed(shard_size, hashes, seed+random::<u64>()))));
        agregated_bf_vector_2.push(Arc::new(Mutex::new(AggregatingBloomFilter::new_with_seed(shard_size, hashes, seed+random::<u64>()))));
    }
    let mut hist_mutex_vector: Vec<Arc<Mutex<u64>>> = Vec::new();
    for _ in 0..nb_files+1{
        hist_mutex_vector.push(Arc::new(Mutex::new(0)));
    }
    reader
        .lines()
        .par_bridge()
        .for_each(|line| {
            let mut bf: BloomFilter = BloomFilter::new_with_seed(size, hashes, seed+1312);
            let filename = line.unwrap();
            println!("{}", filename);
            handle_fasta(filename, &agregated_bf_vector, &agregated_bf_vector_2, &mut bf, modimizer, &hist_mutex_vector);
        });
    Ok(hist_mutex_vector)
}

fn handle_fasta(filename: String, agregated_bf_vector: &Vec<Arc<Mutex<AggregatingBloomFilter>>>, agregated_bf_vector_2: &Vec<Arc<Mutex<AggregatingBloomFilter>>>, bf: &mut BloomFilter, modimizer: u64, hist_mutex_vector: &Vec<Arc<Mutex<u64>>>){
    let mut missing = false;
    let ( reader, _compression) = niffler::get_reader(Box::new(File::open(filename).unwrap())).unwrap();
    let mut fa_reader = Reader::new(reader);
    let mut kmer = RawKmer::<K, KT>::new();
    while let Some(record) = fa_reader.next(){
        let record = record.expect("Error reading record");
        let mut size = 0;
        for (i, nuc) in record.seq().iter().filter_map(KT::from_nuc).enumerate() {
            if size < K - 1{
                kmer = kmer.extend(nuc);
                size += 1;
            }else{
                kmer = kmer.append(nuc);
                let canon = kmer.canonical().to_int();
                if canon%modimizer == 0{
                    missing = bf.insert_if_missing(canon);
                }
                if missing{
                    let mut curr_vec_1 = agregated_bf_vector.get((canon%SHARD_AMOUNT as u64) as usize).unwrap().lock().unwrap();
                    let mut curr_vec_2 = agregated_bf_vector_2.get((canon%SHARD_AMOUNT as u64) as usize).unwrap().lock().unwrap();
                    let count_1 = curr_vec_1.add_and_count(canon/SHARD_AMOUNT as u64);
                    let count_2 = curr_vec_2.add_and_count(canon/SHARD_AMOUNT as u64);
                    drop(curr_vec_1);
                    drop(curr_vec_2);
                    let min_count = min(count_1, count_2);
                    if min_count < hist_mutex_vector.len() as u16{
                        let mut curr_counter = hist_mutex_vector.get(min_count as usize).unwrap().lock().unwrap();
                        *curr_counter += 1;
                        drop(curr_counter);
                        if min_count != 1{
                            let mut prev_counter = hist_mutex_vector.get((min_count-1) as  usize).unwrap().lock().unwrap();
                            *prev_counter -= 1;
                            drop(prev_counter);
                        }
                    }
                    missing = false;
                }
            }
            
        }
    }
}

fn write_output(hist: Vec<Arc<Mutex<u64>>>, nb_files: usize, output: String) -> Result<(), Box<dyn Error>>{

    let mut wtr = Writer::from_path(output)?;
    let header: Vec<u16> = (1..(nb_files+1) as u16).collect();
    let mut hist_to_write = Vec::new();
    wtr.serialize(header)?;
    for i in 0..nb_files+1{
        let val = hist.get(i).unwrap().lock().unwrap();
        hist_to_write.push(val.clone());
        drop(val);
    }
    wtr.serialize(&hist_to_write[1..nb_files+1])?;
    wtr.flush()?;
    Ok(())
}


#[test]
fn test_process_fof_parallel() {
    // Replace "path/to/your/test_file.txt" with the path to your test file
    let filename = "../Data/fof_test.txt";
    let modimizer = 1; // Adjust as needed
    let nb_files = 2; // Adjust as needed
    let size = 1_000_000_000; // Adjust as needed
    let hashes = 1; // Adjust as needed
    let seed = 1515; // Adjust as needed
    let nb_filters = 2; // Adjust as needed
    let shard_amount = 4; // Adjust as needed
    env::set_var("RAYON_NUM_THREADS", "1");
    // Run the function under test
    let result = process_fof_parallel(filename, modimizer, nb_files, size, hashes, seed, nb_filters, shard_amount);

    // Assert that the function returns successfully
    assert!(result.is_ok());

    // Extract the result vector from the mutex
    let hist_mutex = result.unwrap();
    let hist = hist_mutex.lock().unwrap();

    // Add additional assertions based on the expected behavior of your function
    assert_eq!(hist.len(), nb_files + 1);
    assert_eq!(hist[2], 97);
    assert_eq!(hist[0], 0);
    // Add more assertions as needed
}