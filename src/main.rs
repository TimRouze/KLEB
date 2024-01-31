#![allow(dead_code)]

mod bloom;
mod utils;
use clap::Parser;
use seq_io::fasta::Reader;
use std::sync::{Arc, Mutex};
use std::error::Error;
use std::fs::File;
use std::io::{self, BufRead};
use std::cmp::min;//bebou
use std::env;
use::csv::Writer;
use::rayon::prelude::*;
use bloom::{BloomFilter, AggregatingBloomFilter};

const K: usize = 31;
pub type KT = u64;

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Input file (.fasta, .fa)
    input: String,
    /// Output file (defaults to <input>.csv)
    ///#[arg(short, long)]
    ///output: Option<String>,
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
}
fn main() {
    let args = Args::parse();
    let input_fof = args.input.as_str();
    let size =  args.memory * 1_000_000_000;

    env::set_var("RAYON_NUM_THREADS", args.threads.to_string());
    if let Ok(lines_temp) = utils::read_lines(input_fof){
        let nb_files = lines_temp.count();
        match process_fof_parallel(input_fof, args.modimizer, nb_files, size) {
            Ok(hist_mutex) => {
                println!("All {} files have been read...\nWriting output...", nb_files);
                let hist = Arc::try_unwrap(hist_mutex).expect("Failed to unnwrap Arc").into_inner().expect("Failed to get Mutex");
                write_output(hist, nb_files).unwrap();
            }
            Err(err) => eprintln!("Error reading or processing file: {}", err),
        }
    }
//    var |-ma-variable-est-un-kebab-|
//    var ma_variable_est_un_serpent
//    var maVariableEstUnChameaubebou
}

fn process_fof_parallel(filename: &str, modimizer: u64, nb_files: usize, size: usize) -> io::Result<Arc<Mutex<Vec<u64>>>>{
    let file = File::open(filename)?;
    let reader = io::BufReader::new(file);
    let agregated_BF_mutex_1 = Arc::new(Mutex::new(AggregatingBloomFilter::new_with_seed(size, 1, 333)));//bebou
    let agregated_BF_mutex_2 = Arc::new(Mutex::new(AggregatingBloomFilter::new_with_seed(size, 1, 777)));//bebou
    let hist_mutex = Arc::new(Mutex::new(vec![0; nb_files+1]));
    // Process lines in parallel using rayon
    reader
        .lines()
        .par_bridge()
        .for_each(|line| {
            let mut bf: BloomFilter = BloomFilter::new_with_seed(size, 1, 1312);
            let filename = line.unwrap();
            println!("{}", filename);
            handle_fasta(filename, &agregated_BF_mutex_1, &agregated_BF_mutex_2, &mut bf, modimizer, &hist_mutex);
        });

    Ok(hist_mutex)
}

fn handle_fasta(filename: String, agregated_BF_mutex_1: &Arc<Mutex<AggregatingBloomFilter>>, agregated_BF_mutex_2: &Arc<Mutex<AggregatingBloomFilter>>, bf: &mut BloomFilter, modimizer: u64, hist_mutex: &Arc<Mutex<Vec<u64>>>){
    let mut missing = false;
    let ( reader, _compression) = niffler::get_reader(Box::new(File::open(filename).unwrap())).unwrap();
    let mut fa_reader = Reader::new(reader);
    while let Some(record) = fa_reader.next(){
        let record = record.expect("Error reading record");
        for s in record.seq_lines(){
            let seq = String::from_utf8_lossy(s);
            if seq.len() >= 31{
                for _i in 0..(seq.len()-K){
                    let k_mer = utils::str2num(&seq[_i.._i+K]);
                    if k_mer%modimizer == 0{
                        missing = bf.insert_if_missing(utils::canon(k_mer, utils::rev_comp(k_mer)));
                    }
                    if missing{
                        let mut agregated_BF_1 = agregated_BF_mutex_1.lock().unwrap();
                        let mut agregated_BF_2 = agregated_BF_mutex_2.lock().unwrap();
                        let count_1 = agregated_BF_1.add_and_count(utils::canon(k_mer, utils::rev_comp(k_mer)));
                        let count_2 = agregated_BF_2.add_and_count(utils::canon(k_mer, utils::rev_comp(k_mer)));
                        let mut hist = hist_mutex.lock().unwrap();
                        let min_count = min(count_1, count_2);
                        if min_count < hist.len() as u16{
                            hist[min_count as usize] += 1;
                            if min_count != 1{
                                hist[(min_count-1) as usize] -= 1;
                            }
                        }
                        missing = false;
                    }
                }
            }
        }
    }
}

fn write_output(hist: Vec<u64>, nb_files: usize) -> Result<(), Box<dyn Error>>{

    let mut wtr = Writer::from_path("out.csv")?;
    let header: Vec<u16> = (1..(nb_files+1) as u16).collect();
    wtr.serialize(header)?;
    wtr.serialize(&hist[1..(nb_files+1)])?;
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
    env::set_var("RAYON_NUM_THREADS", "1");
    // Run the function under test
    let result = process_fof_parallel(filename, modimizer, nb_files, size);

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