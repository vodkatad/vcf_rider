extern crate vcf_rider;

//use std::io::BufReader;
//use std::fs::File;
use std::env;
use vcf_rider::fasta;
use vcf_rider::pwm;

fn main() {
    let mut args = env::args();
    let binary_name =  args.nth(0).unwrap();
    let fasta_filename = get_next_arg(&mut args, "Missing fasta file argument".to_owned(), &binary_name);
    let pwms_filename = get_next_arg(&mut args, "Missing pwm file argument".to_owned(), &binary_name);
    println!("fasta: {}", fasta_filename);
    println!("pwms: {}", pwms_filename);
    let mut matrixes : Vec<pwm::PWM> = Vec::new();
    if let Ok(pwm_reader) = pwm::PWMReader::open_path(&pwms_filename) {
        for pwm in pwm_reader {
            matrixes.push(pwm);
        }
    }
    for m in matrixes {
        println!("name {}", m.name);
        println!("freq {:?}", m.freq);
    }
    if let Ok(reader) = fasta::FastaReader::open_path(&fasta_filename) {
        for f in reader {
            println!("id: {}", f.id);
            println!("seq: {:?}", f.sequence);
        }
    }
}

fn get_next_arg(args: &mut env::Args, error: String, binary_name: &String) -> String {
    match args.next() {
        Some(x) => x,
        None => {
            println!("Usage: {} <fasta_filename> <pwm_filename>", binary_name);
            println!("{}", error);
            std::process::exit(1);
        }
    }
}
