extern crate vcf_rider;
extern crate rust_htslib;
extern crate itertools;
extern crate bio;

use vcf_rider::rider::*;
use std::env;
use vcf_rider::pwm;
use bio::io::bed;
use std::path::Path;

use rust_htslib::bcf;
use itertools::Itertools;

fn main() {
    let mut args = env::args();
    let binary_name =  args.nth(0).unwrap();
    let vcf_filename = get_next_arg(&mut args, "Missing vcf file argument", &binary_name);
    let pwms_filename = get_next_arg(&mut args, "Missing pwm file argument", &binary_name);
    let bed_filename = get_next_arg(&mut args, "Missing bed file argument", &binary_name);
    let ref_filename = get_next_arg(&mut args, "Missing reference chr argument", &binary_name);
    // TODO argument bg
    //println!("fasta: {}", fasta_filename);
    //println!("pwms: {}", pwms_filename);
    let mut matrixes : Vec<pwm::PWM> = Vec::new();
    if let Ok(pwm_reader) = pwm::PWMReader::open_path(&pwms_filename) {
        for pwm in pwm_reader {
            matrixes.push(pwm);
        }
    }
    let mut pwm_lengths = matrixes.iter().map(|pwm| pwm.get_length());
    // It seems to me that min_max is not there anymore, it is more efficient, if needed the code is in TODO.txt.
    let (min, max) = {
        if let Some(mut min) = pwm_lengths.next() {
            let mut max = min;
            for len in pwm_lengths {
                if len < min {
                    min = len;
                }
                else if len > max {
                    max = len;
                }        
            }
            (min, max)
        } 
        else {
            panic!("No PWM were found in the matrixes file!");
        }        
    };

    // Do not understand: http://hermanradtke.com/2015/06/22/effectively-using-iterators-in-rust.html
    // why here is different?
    // If you want to explore why matrixes.iter.map moves  matrixes to pwm_lengths (its scope) uncomment:
    // for m in matrixes {
    //     println!("{}", m.name);
    // }

    if let Ok(bed_reader) = bed::Reader::from_file(Path::new(&bed_filename)) {
        get_scores(RiderParameters {min_len: min, max_len: max, parameters: &matrixes}, &vcf_filename, bed_reader, &ref_filename);
    } else {
        panic!("Could not open bed file!");
    }

    let vcf_reader = match bcf::Reader::new(&vcf_filename) {
        Ok(x) => x,
        Err(e) => panic!("Could not open vcf file {}", e)
    };
    for sample in vcf_reader.header.samples() {
        let name = sample.into_iter().map(|a| a.to_owned()).map(|x| x as char).collect_vec();
        println!("{:?}", name); // it is ['H', 'G', '0', '0', '1', '0', '5']  -> woooa.
    }
    let mut record = bcf::Record::new();
    loop {
        if let Err(e) = vcf_reader.read(&mut record) {
            if e.is_eof() {
                break;
            } 
        }
        let alleles = record.alleles().into_iter().map(|a| a.to_owned()).collect_vec(); // I do not like this to_owned...
        for allele in alleles[1..].iter() { // why skip the first? The first is the reference. Then are listed all the alternative ones.
            println!("rid2name {:?}" , vcf_reader.header.rid2name(record.rid().unwrap())); // dunno what this is
            println!("record.pos {}" , record.pos() as i32); // Manually checked: the lib converts 1 based coords of vcf to 0 based.
            println!("alleles[0] {}" , alleles[0][0] as char); // this is the reference
            println!("allele {}" , allele[0] as char); // this is the alt allele? Why are they vectors? Because they are encoded as single bases.
            //println!("allele {}" , allele[1] as char); // this will print the second base if the alt allele is for example AC (C). 
        }
        // remember trim_alleles(&mut self)
        
        let genotypes = record.genotypes().unwrap();
        let genotypes = (0..vcf_reader.header.sample_count() as usize).map(|s| {
                        format!("{}", genotypes.get(s))
                        }).collect_vec();
            
        for s in 0..vcf_reader.header.sample_count() as usize {
            println!("genotypes {}", genotypes[s]); // here genotypes[s] is the str (?) 0|0, ok. These are displayable, check their code to understand the inner structure.
            // https://github.com/rust-bio/rust-htslib/blob/master/src/bcf/record.rs
        }
    }
}

// TODO: use a meaningful crate for args management.
fn get_next_arg(args: &mut env::Args, error: &str, binary_name: &str) -> String {
    match args.next() {
        Some(x) => x,
        None => {
            println!("Usage: {} <vcf_filename> <pwm_filename> <bed filename> <fasta ref filename>", binary_name);
            println!("{}", error);
            std::process::exit(1);
        }
    }
}