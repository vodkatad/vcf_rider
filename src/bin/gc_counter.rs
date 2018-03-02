extern crate vcf_rider;

//use std::io::BufReader;
//use std::fs::File;
use std::env;
use vcf_rider::fasta;
use vcf_rider::rider::CanScoreSequence;
use std::collections::hash_map::HashMap;

const C: u8 = 1;
const G: u8 = 2;

pub struct Counter {
    pub nucleos: HashMap<u8, u64>,
    pub name: String
}
impl Counter {
    pub fn new(to_count: Vec<u8>, name: &str) -> Counter {
        let mut nucleos : HashMap<u8, u64> = HashMap::new();
        for n in to_count {
            nucleos.insert(n, 0u64);
        }
        Counter { nucleos: nucleos, name: name.to_owned()}
    }
}

impl CanScoreSequence for Counter {
    fn get_length(&self) -> usize {
        1
    }

    fn get_name(&self) -> &str {
        &self.name
    }

    fn get_score(&self, pos: usize, sequence: &[u8]) -> f64 {
        // we could also keep counts for single nucleotide inside Counter in this case.
        if self.nucleos.contains_key(&sequence[pos]) {
            1f64
        } else {
            0f64
        }
    }
}

fn main() {
    let mut args = env::args();
    let binary_name =  args.nth(0).unwrap();
    let fasta_filename = get_next_arg(&mut args, "Missing fasta file argument".to_owned(), &binary_name);
    let gc_counter = Counter::new(vec!(C, G), &"CG");
    //println!("fasta: {}", fasta_filename);
    //println!("pwms: {}", pwms_filename);
    if let Ok(reader) = fasta::FastaReader::open_path(&fasta_filename) {
        for f in reader {
            //println!("id: {}", f.id);
            //println!("seq: {:?}", f.sequence);
            let mut windows = f.sequence.as_slice().windows(1); 
            // for gc content we need to implement a not sliding but jumping window inside vcf_rider...or jump 1 by 1, will it work? 
            // Super corner case!
            let mut done = false;
            let mut tot = 0f64;
            while !done {
                if let Some(window) = windows.next() {
                    tot += gc_counter.get_score(0usize, window);
                } else {
                    done = true;
                }
            }
            println!("{}\t{}\t{}\t{}", f.id, gc_counter.name, tot, f.sequence.len());
        }
    }
}


fn get_next_arg(args: &mut env::Args, error: String, binary_name: &String) -> String {
    match args.next() {
        Some(x) => x,
        None => {
            println!("Usage: {} <fasta_filename>", binary_name);
            println!("{}", error);
            std::process::exit(1);
        }
    }
}
