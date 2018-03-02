extern crate vcf_rider;
extern crate bio;
extern crate argparse;
//extern crate gc_counter;

use vcf_rider::rider::*;
use bio::io::bed;
use std::path::Path;
use argparse::{ArgumentParser, Store};
//use gc_counter::{G, C, Counter};
//mod gc_counter;
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


// To obtain GC counts for different individuals on regions of interest.
fn main() {
    let mut vcf_filename = "".to_string();
    let mut bed_filename = "".to_string();
    let mut ref_filename = "".to_string();
    let mut associations_filename = "".to_string();
    
    { 
        let mut ap = ArgumentParser::new();
        ap.set_description("Compute TBA on some genomic intervals for the individuals whose mutations are listed in the VCF. Needs a phased vcf. Works on single chromosomes.");
        ap.refer(& mut vcf_filename).add_option(&["-v", "--vcf"], Store, "A vcf for a single chromosome").required();
        ap.refer(& mut bed_filename).add_option(&["-b", "--bed"], Store, "A bed file representing the desired genomic intervals, on a single chromosome").required();
        ap.refer(& mut ref_filename).add_option(&["-r", "--ref"], Store, "A fasta with the reference sequence for the chromosome of interest").required();
        ap.refer(& mut associations_filename).add_option(&["-a", "--assoc"], Store, "Optional filename where bed-mutations associations will be printed.");
        ap.parse_args_or_exit();
    }
    
    let assoc_file_opt = match associations_filename.as_ref() {
        "" => None,
        _ => Some(associations_filename)
    };
    // get_scores is the workhorse of our library, given some RiderParameters (i.e. a set of objects which can "score" a sequence,
    // in this case we want simply to count G and C nucleotides), a vcf filename, a bed reader, a fasta file name with the
    // reference sequence and an optional file where it will print bed entries with their overlapping snps (only for entries with at least one overlapping SNP.)
    // In the end it is not good to print directly results inside the library so it will return an appropriate data structure with results that will be printed here.
    // Right now for our pipelines it is ok to print inside it.
    let gc_counter = Counter::new(vec!(C, G), &"CG");
    let mut vec : Vec<Counter> = Vec::new();
    vec.push(gc_counter);

    if let Ok(bed_reader) = bed::Reader::from_file(Path::new(&bed_filename)) {
        get_scores(RiderParameters {min_len: 1, max_len: 1, parameters: &vec}, &vcf_filename, bed_reader, &ref_filename, assoc_file_opt, false);
    } else {
        panic!("Could not open bed file!");
    }
}