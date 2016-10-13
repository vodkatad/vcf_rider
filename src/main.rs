extern crate vcf_rider;
extern crate bio;

use vcf_rider::rider::*;
use std::env;
use vcf_rider::pwm;
use bio::io::bed;
use std::path::Path;
use std::iter::MinMaxResult::{self, NoElements, OneElement, MinMax};

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
    let pwm_lengths = matrixes.map(|x| x.get_length());
    let min_max = pwm_lengths.min_max();
    rider::get_scores(RiderParameters {matrixes, min_max}, vcf_filename, bed::Reader::from_file(Path::new(bed_filename)));
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