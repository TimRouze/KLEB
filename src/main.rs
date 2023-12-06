#![allow(dead_code)]

mod bloom;
use clap::Parser;
use seq_io::fasta::Reader;
use std::sync::{Arc, Mutex};
use std::error::Error;
use std::fs::File;
use std::path::Path;
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
    #[arg(short = 'M', long)]
    modimizer: bool,
}
fn main() {
    let args = Args::parse();
    let input_fof = args.input.as_str();
    let size =  args.memory * 1_000_000_000;

    env::set_var("RAYON_NUM_THREADS", args.threads.to_string());
    if let Ok(lines_temp) = read_lines(input_fof){
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

fn process_fof_parallel(filename: &str, modimizer: bool, nb_files: usize, size: usize) -> io::Result<Arc<Mutex<Vec<u64>>>>{
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

fn handle_fasta(filename: String, agregated_BF_mutex_1: &Arc<Mutex<AggregatingBloomFilter>>, agregated_BF_mutex_2: &Arc<Mutex<AggregatingBloomFilter>>, bf: &mut BloomFilter, modimizer: bool, hist_mutex: &Arc<Mutex<Vec<u64>>>){
    let mut missing = false;
    let ( reader, _compression) = niffler::get_reader(Box::new(File::open(filename).unwrap())).unwrap();
    let mut fa_reader = Reader::new(reader);
    while let Some(record) = fa_reader.next(){
        let record = record.expect("Error reading record");
        for s in record.seq_lines(){
            let seq = String::from_utf8_lossy(s);
            if seq.len() >= 31{
                for _i in 0..(seq.len()-K){
                    let k_mer = str2num(&seq[_i.._i+K]);
                    if modimizer{
                        if k_mer%2 == 0{
                            missing = bf.insert_if_missing(canon(k_mer, rev_comp(k_mer)));
                        }
                    }else{
                        /* println!("kmer: {}\nrevComp: {}\ncanon: {}", num2str(k_mer), num2str(revcomp), num2str(k_mer_canon));
                        let mut s=String::new();
                        stdin().read_line(&mut s).expect("Did not enter a correct string"); */
                        missing = bf.insert_if_missing(canon(k_mer, rev_comp(k_mer)));
                    }
                    if missing{
                        let mut agregated_BF_1 = agregated_BF_mutex_1.lock().unwrap();
                        let mut agregated_BF_2 = agregated_BF_mutex_2.lock().unwrap();
                        let count_1 = agregated_BF_1.add_and_count(canon(k_mer, rev_comp(k_mer)));
                        let count_2 = agregated_BF_2.add_and_count(canon(k_mer, rev_comp(k_mer)));
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


//bebou
fn canon(k_mer1: u64, k_mer2:u64) -> u64{
    min(k_mer1, k_mer2)
}

fn str2num(k_mer: &str) -> u64{
    let mut res : u64 = 0;
    for character in k_mer.chars(){
        res <<=2;
        res += (character as u64/2)%4;
    }
    res
}

fn num2str(mut k_mer: u64) -> String{
    let mut res = String::from("");
    let mut nuc: u64;
    for _i in 0..K{
        nuc = k_mer%4;
        if nuc == 0{
            res.push('A');
        }else if nuc == 1{
            res.push('C');
        }else if nuc == 2{//bebou
            res.push('T');
        }else if nuc == 3{
            res.push('G');
        }
        k_mer >>= 2;
    }
    res.chars().rev().collect()
}

fn rev_comp(k_mer : u64) -> u64 {
    let mut res = k_mer.reverse_bits();
    res = (res >> 1 & 0x5555_5555_5555_5555) | (res & 0x5555_5555_5555_5555) << 1;
    res ^= 0xAAAA_AAAA_AAAA_AAAA;
    res >> (2 * (32 - K))
}

fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>>
where P: AsRef<Path>, {
    let file = File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}
#[test]
fn test_num2str(){
    let kmer = "AGCTTTTCATTCTGACTGCAACGGGCAATAT";
    let num = str2num(kmer);
    assert_eq!(num2str(num), kmer);
}
#[test]
fn test_rc_64() {
    let kmer = "AGCTTTTCATTCTGACTGCAACGGGCAATAT";
    let revcomp = rev_comp(str2num(kmer));
    let res = num2str(revcomp);
    assert_eq!(res, "ATATTGCCCGTTGCAGTCAGAATGAAAAGCT");
}

#[test]
fn test_canon(){
    let seq = "AGCTTTTCATTCTGACTGCAACGGGCAATAT";
    let kmer = str2num(seq);
    let revcomp = rev_comp(kmer);
    let k_mer_canon = canon(kmer, revcomp);
    assert_eq!(num2str(k_mer_canon), "ATATTGCCCGTTGCAGTCAGAATGAAAAGCT");
}