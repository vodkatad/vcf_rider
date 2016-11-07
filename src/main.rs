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
    let vcf_filename = get_next_arg(&mut args, "Missing vcf file argument", &binary_name);
    let pwms_filename = get_next_arg(&mut args, "Missing pwm file argument", &binary_name);
    let bed_filename = get_next_arg(&mut args, "Missing bed file argument", &binary_name);
    let ref_filename = get_next_arg(&mut args, "Missing reference chr argument", &binary_name);
    let bg = vec!(0.298947240099661, 0.200854143743417, 0.200941012710477, 0.299257603446445);
    // TODO argument bg
    //println!("fasta: {}", fasta_filename);
    //println!("pwms: {}", pwms_filename);
    let mut matrixes : Vec<pwm::PWM> = Vec::new();
    if let Ok(pwm_reader) = pwm::PWMReader::open_path(&pwms_filename) {
        for mut pwm in pwm_reader {
            pwm.compute_ll(&bg);
            matrixes.push(pwm);
        }
    }
    // pwm_lengths is an Iter, matrixes is borrowed until end of scope. Use collect() to avoid.
    // but without collect and with .next() instead of pop() it's more efficient?
    let mut pwm_lengths : Vec<usize> = matrixes.iter().map(|pwm| pwm.get_length()).collect();
    // It seems to me that min_max is not there anymore, it is more efficient, if needed the code is in TODO.txt.
    let (min, max) = {
        if let Some(mut min) = pwm_lengths.pop() {
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
    
    if let Ok(bed_reader) = bed::Reader::from_file(Path::new(&bed_filename)) {
        get_scores(RiderParameters {min_len: min, max_len: max, parameters: &matrixes}, &vcf_filename, bed_reader, &ref_filename);
    } else {
        panic!("Could not open bed file!");
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