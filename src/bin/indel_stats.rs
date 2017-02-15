extern crate vcf_rider;
extern crate bio;

use vcf_rider::rider::*;
use std::env;
use bio::io::bed;
use std::path::Path;
//use vcf_rider::fasta;
use vcf_rider::mutations;
use std::collections::VecDeque;


fn main() {
    let mut args = env::args();
    let binary_name =  args.nth(0).unwrap();
    let vcf_filename = get_next_arg(&mut args, "Missing vcf file argument", &binary_name);
    let bed_filename = get_next_arg(&mut args, "Missing bed file argument", &binary_name);
    if let Ok(mut bed_reader) = bed::Reader::from_file(Path::new(&bed_filename)) {
        if let Ok(vcf) = mutations::VcfReader::open_path(&vcf_filename) {
            let mut vcf_reader = vcf;
            /*for sample in & vcf_reader.samples {
                println!("sample {}", sample);
            }
            for snp in vcf {
                println!("snp {:?} {:?} {:?} {} {:?}", snp.pos, snp.sequence_ref, snp.sequence_alt, snp.id, snp.genotypes);
            }*/
            
            let n_samples = vcf_reader.samples.len();
            let mut snps_buffer : VecDeque<mutations::Mutation> = VecDeque::new();
            for r in bed_reader.records() {
                let record = r.ok().expect("Error reading record");
                //println!("bed name: {}", record.name().expect("Error reading name"));
                // chr check
                // add real windows (length 30) and compute n. of samples in different groups foreach position.
                let window = mutations::Coordinate{chr: "".to_owned(), start: record.start(), end: record.end()};
                let n_overlapping = find_overlapping_snps(& window, &mut vcf_reader, &mut snps_buffer);          
                if snps_buffer.iter().any(|x| x.is_indel) { // filter_map
                    let n_indel = snps_buffer.iter().fold(0, |sum, x| if x.is_indel { sum + 1 } else { sum } );
                    println!("{}\t{}\t{}", record.name().expect("Error reading name"), n_overlapping, n_indel);
                }
            }
        }
    } else {
        panic!("Could not open bed file!");
    }
}

// TODO: use a meaningful crate for args management.
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