extern crate vcf_rider;

use std::env;
use vcf_rider::mutations;

fn main() {
    let mut args = env::args();
    let binary_name =  args.nth(0).unwrap();
    let vcf_filename = get_next_arg(&mut args, "Missing vcf file argument", &binary_name);
    if let Ok(vcf) = mutations::VcfReader::open_path(&vcf_filename, true) {
        for sample in &vcf.samples {
            println!("sample {}", sample);
        }
        for snp in vcf {
            println!("snp {:?} {:?} {:?} {} {:?}", snp.pos, snp.sequence_ref, snp.sequence_alt, snp.id, snp.genotypes);
        }
        
    } else {
        panic!("Could not open vcf file!");
    }
}

fn get_next_arg(args: &mut env::Args, error: &str, binary_name: &str) -> String {
    match args.next() {
        Some(x) => x,
        None => {
            println!("Usage: {} <vcf_filename> <bed filename>", binary_name);
            println!("{}", error);
            std::process::exit(1);
        }
    }
}