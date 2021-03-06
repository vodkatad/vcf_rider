extern crate vcf_rider;

use std::env;
use vcf_rider::fasta;
use vcf_rider::rider::CanScoreSequence;
use vcf_rider::pwm;

// This binary is used to compute TBA on a fasta file, it was just used as a first example
// to compare the speed of this rust implementation with matrix_rider (C).
// They proved to be very similar.
// Warning: this uses a fixed bg defined in fasta.rs.
fn main() {
    let mut args = env::args();
    let binary_name =  args.nth(0).unwrap();
    let fasta_filename = get_next_arg(&mut args, "Missing fasta file argument".to_owned(), &binary_name);
    let pwms_filename = get_next_arg(&mut args, "Missing pwm file argument".to_owned(), &binary_name);
    //println!("fasta: {}", fasta_filename);
    //println!("pwms: {}", pwms_filename);
    let mut matrixes : Vec<pwm::PWM> = Vec::new();
    if let Ok(pwm_reader) = pwm::PWMReader::open_path(&pwms_filename) {
        for pwm in pwm_reader {
            matrixes.push(pwm);
        }
    }
    if let Ok(reader) = fasta::FastaReader::open_path(&fasta_filename) {
        for f in reader {
            //println!("id: {}", f.id);
            //println!("seq: {:?}", f.sequence);

            for i in 0..matrixes.len() { // the iterator here does not work
                let mut m = matrixes.get_mut(i).unwrap();
                //println!("name {}", m.name);
                //println!("freq {:?}", m.freq);
                //println!("len {}", m.freq.len());
                m.compute_ll(&f.background);
                //println!("freq {:?}", m.ll);
                //println!("freq {:?}", m.llrc);
                let mut windows = f.sequence.as_slice().windows(m.get_length());
                let mut done = false;
                let mut tot = 0f64;
                while !done {
                    if let Some(window) = windows.next() {
                        tot += m.get_score(0usize, window);
                    } else {
                        done = true;
                    }
                }
                println!("{}\t{}\t{}", f.id, m.name, tot);
            }
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
