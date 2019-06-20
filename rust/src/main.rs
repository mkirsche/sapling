use seq_io::fasta::{self};
use failure::Error;
use std::env;
use std::io::Write;
use std::fs::File;
use std::str;
use suffix::SuffixTable;

fn read_genome(filename : &str) -> Result<String, Error> {
    let mut reader = fasta::Reader::from_path(filename)?;
    let mut total_seq = "".to_string();
    while let Some(record) = reader.next() {
        let owned = record?.to_owned_record();
        let cur_seq : &str = &str::from_utf8(&owned.seq)?;
        for c in cur_seq.chars() {
            if c == 'A' || c == 'C' || c == 'G' || c == 'T' {
                total_seq.push(c);
            }
        }
    }
    Ok(total_seq)
}

fn suffix_array(s : &str) -> Vec<u32> {
    let st = SuffixTable::new(s);
    let res_array = st.table();
    res_array.to_vec()
}

fn main() {
    let args : Vec<String> = env::args().collect();
    if args.len() != 3 {
        println!("Usage: cargo run <fastafile> <outputfile>");
        return;
    } 
    let filename : &str = &args[1].as_str();
    let outfilename : &str = &args[2].as_str();
    let k : usize = 21;
    let genome = read_genome(filename).unwrap();
    println!("Read genome of length {:?}", genome.len());
    let sa = suffix_array(&genome);
    println!("Constructed suffix array");
    let mut writer = File::create(outfilename).unwrap();
    for (_, &val) in sa.iter().enumerate() {
        let v = val as usize;
        let substring : String = (&genome)[v..(v + k).min(genome.len())].to_string();
        writeln!(&mut writer, "{:?} {}", val, substring).unwrap();
    }
}
