extern crate vcf_rider;
extern crate bio;
extern crate argparse;

use vcf_rider::rider::*;
use vcf_rider::pwm;
use bio::io::bed;
use std::path::Path;
use argparse::{ArgumentParser, Store};
use std::io::BufReader;
use std::fs::File;
use std::io::BufRead;

// TODO move this to bin? This could be called vcf_rider_tba
fn main() {
    let mut vcf_filename = "".to_string();
    let mut pwms_filename = "".to_string(); 
    let mut bed_filename = "".to_string();
    let mut ref_filename = "".to_string();
    let mut associations_filename = "".to_string();
    let mut bg_filename = "".to_string();
    { 
        let mut ap = ArgumentParser::new();
        ap.set_description("Compute TBA on some genomic intervals for the individuals whose mutations are listed in the VCF. Needs a phased vcf. Works on single chromosomes.");
        ap.refer(& mut vcf_filename).add_option(&["-v", "--vcf"], Store, "A phased vcf for a single chromome").required();
        ap.refer(& mut pwms_filename).add_option(&["-p", "--pwm"], Store, "PWM file in the format required by matrix rider (name, pos, a, c, g, t), with counts, no zeroes.").required();
        ap.refer(& mut bed_filename).add_option(&["-b", "--bed"], Store, "A bed file representing the desired genomic intervals, on a single chromosome").required();
        ap.refer(& mut ref_filename).add_option(&["-r", "--ref"], Store, "A fasta with the reference sequence for the chromosome of interest").required();
        ap.refer(& mut bg_filename).add_option(&["-f", "--bg"], Store, "A tab delimited single line (a, c, g, t) with the background frequencies to be used").required();
        ap.refer(& mut associations_filename).add_option(&["-a", "--assoc"], Store, "Optional filename where bed-mutations associations will be printed.");
        ap.parse_args_or_exit();
    }
    //println!("fasta: {}", ref_filename);
    //println!("pwms: {}", pwms_filename);
    let mut read_bg = match File::open(&bg_filename) {
        Ok(file) => BufReader::new(file),
        Err(_) =>  panic!("Error while reading bg file, does it have a single tab separated line?")
    };
    let mut bg = String::new();
    let bg_vec = match read_bg.read_line(& mut bg) {
        Ok(_) => bg.trim_right().split("\t").to_owned(),
        Err(_) => panic!("Error while reading bg file, does it have a single tab separated line?")
    };
    let float_bg = bg_vec.map(|x| { x.parse::<f64>().unwrap()}).collect();
    let mut matrixes : Vec<pwm::PWM> = Vec::new();
    if let Ok(pwm_reader) = pwm::PWMReader::open_path(&pwms_filename) {
        for mut pwm in pwm_reader {
            pwm.compute_ll(&float_bg);
            matrixes.push(pwm);
        }
    } else {
        panic!("Error while reading PWM file!")
    }

    let assoc_file_opt = match associations_filename.as_ref() {
        "" => None,
        _ => Some(associations_filename)
    };
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
    // get_scores is the workhorse of our library, given some RiderParameters (i.e. a set of objects which can "score" a sequence,
    // in this case PWMs, and their minimum and maximum needed sequence lengths), a vcf filename, a bed reader, a fasta file name with the
    // reference sequence and an optional file where it will print bed entries with their overlapping snps (only for entries with at least one overlapping SNP.)
    // In the end it is not good to print directly results inside the library so it will return an appropriate data structure with results that will be printed here.
    // Right now for our pipelines it is ok to print inside it.
    if let Ok(bed_reader) = bed::Reader::from_file(Path::new(&bed_filename)) {
        get_scores(RiderParameters {min_len: min, max_len: max, parameters: &matrixes}, &vcf_filename, bed_reader, &ref_filename, assoc_file_opt);
    } else {
        panic!("Could not open bed file!");
    }
}