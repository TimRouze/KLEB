#![allow(dead_code)]

mod bloom;
use clap::Parser;
use seq_io::fasta::{Reader};
use std::os::linux::net::SocketAddrExt;
use std::thread;
use std::sync::{Arc, Mutex};
use std::error::Error;
use std::fs::File;
use std::path::Path;
use std::io::{self, BufRead, stdin};
use std::cmp::min;//bebou
use::csv::Writer;
use niffler;
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
    #[arg(short, long)]
    threads: Option<usize>,
    /// Memory (in GB) allocated to Bloom filters (defaults to input size)
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
    let modimizer = if args.modimizer {
        true
    } else {
        false
    };
    let threads = if let Some(num) = args.threads {
        num
    } else {
        std::thread::available_parallelism().map_or(1, |x| x.get())
    };
    /*let output_filename = if let Some(filename) = args.output{
        filename
    }else if let Some((begin, end)) = input_filename.rsplit_once('.'){
        begin.to_owned() + ".csv" + end
    } else{
        input_filename.to_owned() + ".csv"
    };
    */
    let size =  args.memory * 1_000_000_000 / 2;
    if let Ok(lines_temp) = read_lines(input_fof){
        let nb_files = lines_temp.count();
        println!("{}", nb_files);
        let mut aBF_mutex = Arc::new(Mutex::new(AggregatingBloomFilter::new_with_seed(size, args.hashes, args.seed)));//bebou
        let mut hist_mutex = Arc::new(Mutex::new(vec![0; (nb_files+1)]));
        if let Ok(lines) = read_lines(input_fof){
            let mut threads = vec![];
            for line in lines{
                println!("{:?}",line);
                let aBF_mutex_clone = Arc::clone(&aBF_mutex);
                let hist_mutex_clone = Arc::clone(&hist_mutex);
                threads.push(thread::spawn(move|| {
                    if let Ok(filename) = line{
                        let mut aBF_mut = aBF_mutex_clone.lock().unwrap();
                        let mut hist_mut = hist_mutex_clone.lock().unwrap();
                        let mut bf: BloomFilter = BloomFilter::new_with_seed(size, args.hashes, args.seed);
                        handle_fasta(filename, &mut aBF_mut, &mut bf, modimizer, &mut hist_mut);
                    }
                }));
            }
            for handle in threads {
                let _ = handle.join();
            }
            println!("All {} files have been read...\nWriting output...", nb_files);
            //aBF.clear();
            let hist = Arc::try_unwrap(hist_mutex).unwrap();
            write_output(hist.into_inner().unwrap(), nb_files).unwrap();
            
        }
    }
//    var |-ma-variable-est-un-kebab-|
//    var ma_variable_est_un_serpent
//    var maVariableEstUnChameaubebou
}

fn handle_fasta(filename: String, aBF: &mut AggregatingBloomFilter, bf: &mut BloomFilter, modimizer: bool, hist: &mut Vec<u16>){
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
                        let count = aBF.add_and_count(canon(k_mer, rev_comp(k_mer)));
                        if count < hist.len() as u16{
                            hist[count as usize] += 1;
                            hist[(count-1) as usize] -= 1;
                        }
                        missing = false;
                    }
                }
            }
        }
    }
}

fn write_output(hist: Vec<u16>, nb_files: usize) -> Result<(), Box<dyn Error>>{
    /*let mut result_vec: Vec<u64> = vec![0; nb_files as usize];
    let mut counter = 0;
    let mut erroneous_count = 0;
    for elem in aBF.counts{
        if elem != 0 && elem <= nb_files{
            counter += 1;
            result_vec[(elem-1) as usize] += 1;
            /*println!("Elem = {}\nvec[elem] = {}\ncounter = {}", elem-1, result_vec[(elem-1) as usize], counter);
            let mut s=String::new();
            stdin().read_line(&mut s).expect("Did not enter a correct string");*/
        }else if elem > nb_files{
            erroneous_count += 1;
        }
    }
    println!("I have seen: {} k_mers", counter);
    println!("I have seen: {} collisions", erroneous_count);
    println!("Total number of k_mers seen is: {}", (counter+erroneous_count));*/
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
    return res;
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