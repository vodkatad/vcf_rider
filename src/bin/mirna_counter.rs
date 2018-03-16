extern crate vcf_rider;

use std::env;
use vcf_rider::fasta;
use vcf_rider::rider::CanScoreSequence;
use vcf_rider::mirna;

// This binary is used to compute TBA on a fasta file, it was just used as a first example
// to compare the speed of this rust implementation with matrix_rider (C).
// They proved to be very similar.
// Warning: this uses a fixed bg defined in fasta.rs.
fn main() {
    let mut args = env::args();
    let binary_name =  args.nth(0).unwrap();
    let fasta_filename = get_next_arg(&mut args, "Missing fasta file argument".to_owned(), &binary_name);
    let mirna_filename = get_next_arg(&mut args, "Missing mirna file argument".to_owned(), &binary_name);
    
    let mut mirna : Vec<mirna::Seed> = Vec::new();
    if let Ok(seed_reader) = mirna::SeedReader::open_path(&mirna_filename) {
        for mir in seed_reader {
            mirna.push(mir);
            println!("{}, {:?}", mir.name, mir.sequence);
        }
    }
/*    if let Ok(reader) = fasta::FastaReader::open_path(&fasta_filename) {
        for f in reader {
            //println!("id: {}", f.id);
            //println!("seq: {:?}", f.sequence);

                println!("{}\t{}\t{}", f.id, m.name, tot);
            }
        }

    }*/
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
