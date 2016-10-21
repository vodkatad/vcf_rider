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
    //  m.compute_ll(&f.background); -> has to be done before calling get_scores

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