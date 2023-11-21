#![allow(dead_code)]

mod bloom;
use clap::Parser;
use seq_io::fasta::{Reader};
use std::error::Error;
use std::fs::{metadata, File};
use std::path::Path;
use std::io::{self, BufRead, stdin};
use std::cmp::min;
use std::collections::HashMap;//bebou
use::csv::Writer;
use bloom::{BloomFilter, AggregatingBloomFilter};


const K: u8 = 31;

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Input file (.fasta, .fa)
    input: String,
    /// Output file (defaults to <input>.csv)
    ///#[arg(short, long)]
    ///output: Option<String>,
    /// Memory (in MB) allocated to Bloom filters (defaults to input size)
    #[arg(short, long, default_value_t = 4)]
    memory: usize,
    /// Number of hashes used in Bloom filters
    #[arg(short = 'H', long, default_value_t = 3)]
    hashes: usize,
    /// Seed used for hash functions
    #[arg(short, long, default_value_t = 101010)]
    seed: u64,
}
fn main() {
    let args = Args::parse();
    let input_fof = args.input.as_str();
    /*let output_filename = if let Some(filename) = args.output{
        filename
    }else if let Some((begin, end)) = input_filename.rsplit_once('.'){
        begin.to_owned() + ".csv" + end
    } else{
        input_filename.to_owned() + ".csv"
    };
    */
    //let mut nb_files : u32 = 0;
    let size =  args.memory * 1_000_000 / 2;
    let mut bf: BloomFilter = BloomFilter::new_with_seed(size, args.hashes, args.seed);                
    let mut aBF: AggregatingBloomFilter = AggregatingBloomFilter::new_with_seed(size, args.hashes, args.seed);//bebou
    if let Ok(lines) = read_lines(input_fof){
        for line in lines{
            if let Ok(filename) = line{
                //nb_files += 1;
                println!("{}", filename);
                let mut reader = Reader::from_path(&filename).unwrap();
                while let Some(record) = reader.next(){
                    let record = record.expect("Error reading record");
                    for s in record.seq_lines(){
                        let seq = String::from_utf8_lossy(s);
                        for i in 0..(seq.len()-K as usize){
                            let k_mer = str2num(&seq[i..i+K as usize]);
                            let revcomp = rev_comp(k_mer);
                            let k_mer_canon = canon(k_mer, revcomp);
                            /* println!("kmer: {}\nrevComp: {}\ncanon: {}", num2str(k_mer), num2str(revcomp), num2str(k_mer_canon));
                            let mut s=String::new();
                            stdin().read_line(&mut s).expect("Did not enter a correct string"); */
                            let missing = bf.insert_if_missing(k_mer_canon);
                            if missing{
                                //TODO CHECK EXACTLY WHAT ADD DOES
                                aBF.add(k_mer_canon);
                            }
                            /*
                            en récupérant la valeur de retour de "insert_if_missing"
                            On pourrait faire l'aggregation ici, a la volée:
                            if missing
                                aBF.insert_or_increment(k_mer_canon)
                             */
                        }
                        bf.clear();
                    }
                }
            }
        }
        println!("All files have been read...\nWriting output...");
        //aBF.clear();
        write(aBF.return_non_zero()).unwrap();
    }
//    var |-ma-variable-est-un-kebab-|
//    var ma_variable_est_un_serpent
//    var maVariableEstUnChameaubebou
}

fn write(non_zeros: Vec<u16>) -> Result<(), Box<dyn Error>>{
    let mut results = HashMap::new();
    for i in 0..non_zeros.len(){
        let val = results.entry(non_zeros[i]).or_insert(1);
        *val += 1;
    }
    let mut wtr = Writer::from_writer(io::stdout());
    for elem in results{
        wtr.serialize((elem.0, elem.1))?;
    }
    
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
    for i in 0..K{
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