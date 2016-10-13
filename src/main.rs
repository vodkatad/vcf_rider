extern crate vcf_rider;
extern crate bio;

use vcf_rider::rider::*;
use std::env;
use vcf_rider::pwm;
use bio::io::bed;
use std::path::Path;

fn main() {
    let mut args = env::args();
    let binary_name =  args.nth(0).unwrap();
    let vcf_filename = get_next_arg(&mut args, "Missing vcf file argument".to_owned(), &binary_name);
    let pwms_filename = get_next_arg(&mut args, "Missing pwm file argument".to_owned(), &binary_name);
    let bed_filename = get_next_arg(&mut args, "Missing bed file argument".to_owned(), &binary_name);
    //println!("fasta: {}", fasta_filename);
    //println!("pwms: {}", pwms_filename);
    let mut matrixes : Vec<pwm::PWM> = Vec::new();
    if let Ok(pwm_reader) = pwm::PWMReader::open_path(&pwms_filename) {
        for pwm in pwm_reader {
            matrixes.push(pwm);
        }
    }
    //let mut pwm_lengths = matrixes.map(|x| x.get_length());
    let mut min = usize::max_value();
    let mut max = 0usize;
    for i in 0..matrixes.len() { 
        let len = matrixes.get(i).unwrap().get_length();
        if len < min {
            min = len;
        }
        else if len > max {
            max = len
        }        
    }
    if let Ok(bed_reader) = bed::Reader::from_file(Path::new(&bed_filename)) {
        get_scores(RiderParameters {min_len: min, max_len: max, parameters: matrixes}, &vcf_filename, bed_reader);
    } else {
        panic!("Could not open bed file!");
    }
}

// TODO: use a meaningful crate for args management.
fn get_next_arg(args: &mut env::Args, error: String, binary_name: &String) -> String {
    match args.next() {
        Some(x) => x,
        None => {
            println!("Usage: {} <vcf_filename> <pwm_filename> <bed filename>", binary_name);
            println!("{}", error);
            std::process::exit(1);
        }
    }
}