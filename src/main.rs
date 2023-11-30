#![allow(dead_code)]

mod bloom;
//mod dashbloom;
//mod kmer;
//mod reads;
use clap::Parser;
use seq_io::fasta::{Reader};
use std::error::Error;
use std::fs::File;
use std::path::Path;
use std::io::{self, BufRead, stdin};
use std::cmp::min;//bebou
use::csv::Writer;
use niffler;
use bloom::{BloomFilter, AggregatingBloomFilter};
//use reads::{BaseRecord, Fasta, ReadProcess};
//use dashbloom::{BloomFilter, AggregatingBloomFilter};
//use kmer::{Base, Kmer, RawKmer};



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
    //env::set_var("RUST_BACKTRACE", "full");
    let args = Args::parse();
    let input_fof = args.input.as_str();
    let modimizer = if args.modimizer {
        true
    } else {
        false
    };
    /*let output_filename = if let Some(filename) = args.output{
        filename
    }else if let Some((begin, end)) = input_filename.rsplit_once('.'){
        begin.to_owned() + ".csv" + end
    } else{
        input_filename.to_owned() + ".csv"
    };
    */
    //let mut nb_files : u32 = 0;
    let size =  args.memory * 1_000_000_000;
    let mut bf: BloomFilter = BloomFilter::new_with_seed(size, args.hashes, args.seed);                
    let mut aBF: AggregatingBloomFilter = AggregatingBloomFilter::new_with_seed(size, args.hashes, args.seed);//bebou
    let mut counter = 0;
    let mut nb_files = 0;
    let mut missing = false;
    if let Ok(lines) = read_lines(input_fof){
        for line in lines{
            if let Ok(filename) = line{
                println!("{}",filename);
                nb_files += 1;
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
                                        if missing{
                                            aBF.add(canon(k_mer, rev_comp(k_mer)));
                                        }
                                    }
                                }else{
                                    /* println!("kmer: {}\nrevComp: {}\ncanon: {}", num2str(k_mer), num2str(revcomp), num2str(k_mer_canon));
                                    let mut s=String::new();
                                    stdin().read_line(&mut s).expect("Did not enter a correct string"); */
                                    //bf.insert(canon(k_mer, rev_comp(k_mer)));
                                    missing = bf.insert_if_missing(canon(k_mer, rev_comp(k_mer)));
                                    if missing{
                                        aBF.add(canon(k_mer, rev_comp(k_mer)));
                                    }
                                }
                                if missing{
                                    counter += 1;
                                    missing = false;
                                }
                            }
                        }
                    }
                }
                bf.clear();
            }
        }
        println!("All {} files have been read...\nWriting output...", nb_files);
        //aBF.clear();
        write(aBF, nb_files).unwrap();
        println!("nb unique k_mers: {}", counter/nb_files);
    }
//    var |-ma-variable-est-un-kebab-|
//    var ma_variable_est_un_serpent
//    var maVariableEstUnChameaubebou
}


fn write(aBF: AggregatingBloomFilter, nb_files: u16) -> Result<(), Box<dyn Error>>{
    let mut result_vec: Vec<u64> = vec![0; nb_files as usize];
    let mut counter = 0;
    let mut erroneous_count = 0;
    for elem in aBF.counts{
        if elem != 0 && elem <= nb_files{
            counter += 1;
            result_vec[(elem-1) as usize] += 1;
            /*println!("Elem = {}\nvec[elem] = {}\ncounter = {}", elem, result_vec[(elem) as usize], counter);
            let mut s=String::new();
            stdin().read_line(&mut s).expect("Did not enter a correct string");*/
        }else if elem > nb_files{
            erroneous_count += 1;
        }
    }
    println!("I have seen: {} k_mers", counter);
    println!("I have seen: {} collisions", erroneous_count);
    println!("Total number of k_mers seen is: {}", (counter+erroneous_count));
    let mut wtr = Writer::from_path("out.csv")?;
    let header: Vec<u16> = (1..(nb_files+1)).collect();
    wtr.serialize(header)?;
    wtr.serialize(result_vec)?;
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