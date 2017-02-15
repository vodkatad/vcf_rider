extern crate vcf_rider;
extern crate bio;
extern crate itertools; 

use vcf_rider::rider::*;
use std::env;
use bio::io::bed;
use std::path::Path;
//use vcf_rider::fasta;
use vcf_rider::mutations;
use std::collections::VecDeque;
use itertools::Itertools;

fn main() {
    let mut args = env::args();
    let binary_name =  args.nth(0).unwrap();
    let vcf_filename = get_next_arg(&mut args, "Missing vcf file argument", &binary_name);
    let bed_filename = get_next_arg(&mut args, "Missing bed file argument", &binary_name);
    if let Ok(mut bed_reader) = bed::Reader::from_file(Path::new(&bed_filename)) {
        if let Ok(vcf) = mutations::VcfReader::open_path(&vcf_filename, true) {
            let mut vcf_reader = vcf;
            /*for sample in & vcf_reader.samples {
                println!("sample {}", sample);
            }
            for snp in vcf {
                println!("snp {:?} {:?} {:?} {} {:?}", snp.pos, snp.sequence_ref, snp.sequence_alt, snp.id, snp.genotypes);
            }*/
            
            let n_samples = vcf_reader.samples.len();
            let mut snps_buffer : VecDeque<mutations::Mutation> = VecDeque::new();
            // We initialize a vector that as indexes has chr_samples ids (chr_M sample 1: 1, chr_P sample 1: 2 ..) and
            // as values the group id.
            for r in bed_reader.records() {
                let mut groups : Vec<u32> =  vec![0; n_samples*2]; // Vec::with_capacity(n_samples*2);
                let record = r.ok().expect("Error reading record");
                //println!("bed name: {}", record.name().expect("Error reading name"));
                // chr check
                // add real windows (length 30) and compute n. of samples in different groups foreach position.
                let window = mutations::Coordinate{chr: "".to_owned(), start: record.start(), end: record.end()};
                let n_overlapping = find_overlapping_snps(& window, &mut vcf_reader, &mut snps_buffer);      
                let n_indel = count_groups(&snps_buffer, n_overlapping, &mut groups, n_samples);
                let n_groups = groups.iter().max().unwrap()+1;
                let group_structure = groups.into_iter().enumerate().map(|x| 
                                                        format!("{:?}:{}", x.0, x.1)).join(";"); // sample : group ; sample:group.
                println!("{}\t{}\t{}\t{}\t{}", record.name().expect("Error reading name"), n_overlapping, n_indel, n_groups, group_structure)
            }
        }
    } else {
        panic!("Could not open bed file!");
    }
}

fn count_groups(snps_buffer: & VecDeque<mutations::Mutation>, n_overlapping: u32, groups: &mut Vec<u32>, n_samples: usize) -> u32 {
    let mut n_indel = 0;    
    for (i_snp, snp) in snps_buffer.iter().enumerate() {
        if snp.is_indel && i_snp <= n_overlapping as usize { // i > n_overlapping we have finished the overlapping snps (the last one is just waiting in the buffer)
            n_indel += 1;
            if snp.genotypes.iter().any(|x| x.0 || x.1) {
                // we have a bisection
                //n_groups = n_groups * 2;
                // even groups have no indels at this run.
                for i_sample in 0 .. n_samples {
                    let allele = snp.genotypes[i_sample];
                    let old_group_0 = groups[i_sample];
                    let old_group_1 = groups[i_sample+n_samples];
                    groups[i_sample] = match allele.0 {
                        true => match old_group_0 {
                            0 => 1,
                            _ => 2u32.pow(old_group_0)+1
                        },
                        false => match old_group_0 {
                            0 => 0,
                            _ => 2u32.pow(old_group_0)
                        }
                    };
                    groups[i_sample+n_samples] = match allele.1 {
                        true => match old_group_1 {
                            0 => 1,
                            _ => 2u32.pow(old_group_1)+1,
                        },
                        false => match old_group_1 {
                            0 => 0,
                            _ => 2u32.pow(old_group_1)
                        }
                    };
                    /*if allele.0 {
                        groups[i_sample] = 2u32.pow(old_group_0)+1;
                    } else {
                        groups[i_sample] = 2u32.pow(old_group_0);
                    }
                    if allele.1 { // ouf, bad, debug then refactor
                        groups[i_sample+n_samples] = 2u32.pow(old_group_1)+1;
                    } else {
                        groups[i_sample+n_samples] = 2u32.pow(old_group_1);
                    }
                    */
                }
            }
        }
    }
    n_indel
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