use rust_htslib::bcf;
use std::io;
use itertools::Itertools;

#[derive(Debug)]
pub struct Coordinate {
    pub chr: String,
    pub start: u64, // 0 based, end inclusive
    pub end: u64
}

// use cmp and Ordering::Less? 
#[derive(Eq, PartialEq, Debug)]
pub enum Position {
    Before,
    Overlapping,
    After
}
impl Coordinate {
    /// Returns the position of self relative to other.
    // right now not using chr
    pub fn relative_position(&self, other : &Coordinate) -> Position {
         if (self.start >= other.start && self.start < other.end) ||
                (self.end > other.start && self.end <= other.end) ||
                (self.start < other.start && self.end >= other.end) {
            Position::Overlapping
        }
        else if self.end <= other.start {
            Position::Before
        } else {
            Position::After
        }
                
    }
}

#[derive(Debug)]
pub struct Mutation {
    pub id: String,
    pub pos: Coordinate,
    pub sequence_ref: Vec<u8>,
    pub sequence_alt: Vec<u8>,
    pub genotypes: Vec<(bool, bool)>
}

pub struct VcfReader {
    reader: bcf::Reader,
    pub samples: Vec<String>
}

// Adding a layer of abstration, I am not sure that we will use the lib.
impl VcfReader {
    pub fn open_path(path: &str) -> io::Result<VcfReader> {
        match bcf::Reader::new(&path) {
            Ok(reader) => {
                let samples = reader.header.samples().into_iter().map(|sample| String::from_utf8(sample.to_owned()).unwrap()).collect();
                Ok(VcfReader { reader: reader, samples: samples })
            }
            Err(_) => Err(io::Error::last_os_error())
            // How do errors work? rust_htslib::bcf::BCFError
        }
    }
}

fn get_sequence(seq : &Vec<u8>) -> Vec<u8> {
    let mut res : Vec<u8> = Vec::<u8>::with_capacity(1);
    for nuc in seq {
        res.push(match *nuc {
            b'A' => 0u8,
            b'C' => 1u8,
            b'G' => 2u8,
            b'T' => 3u8,
            b'N' => 4u8,
            _ => panic!("Vcf with a not allowed char {}", *nuc as char),
        });
    }
    if res.len() > 1 {
        panic!("Right now we do not handle more than SNPs!");
    }
    res
}
fn decode_allele(c: char) -> bool {
    match c {
        '0' => false,
        '1' => true,
         x => panic!("Wrongly encoded snp in vcf {:?}", x)
    }
}

fn decode_genotype(geno: String) -> (bool, bool) {
    // this will need a lot of work inside the lib
    // here genotypes[s] is the str (?) 0|0, ok. These are displayable, check their code to understand the inner structure.
    // https://github.com/rust-bio/rust-htslib/blob/master/src/bcf/record.rs
    let mut genos = geno.chars();
    let a1 = genos.next().unwrap();
    if a1 == '.' {
        return (true, true) // missing data as reference.
    }
    let sep = genos.next().unwrap();
    let a2 = genos.next().unwrap();
    if sep != '|' {
        panic!("Only phased genotypes are supported! {}", sep)
    }
    (decode_allele(a1),decode_allele(a2))
}

impl Iterator for VcfReader {
    type Item = Mutation;

    fn next(&mut self) -> Option<Mutation> {
        let mut record : bcf::Record = bcf::Record::new();
        if let Ok(_) = self.reader.read(&mut record) {
            let id = "".to_owned();
            let chr = "?".to_owned(); // where is it? rid?
            let coord = record.pos() as u64; // Manually checked: the lib converts 1 based coords of vcf to 0 based. bed records have u64
            let alleles = record.alleles().into_iter().map(|a| a.to_owned()).collect_vec(); // I do not like this to_owned...
            let refe : Vec<u8> = get_sequence(&alleles[0]);
            let mut alt  : Vec<u8> = Vec::<u8>::with_capacity(1);
            let mut found_alt = 0;
            for allele in alleles[1..].iter() { // why skip the first? The first is the reference. Then are listed all the alternative ones.
                //println!("allele {}" , allele[0] as char); // this is the alt allele? Why are they vectors? Because they are encoded as single bases.
                alt.push(allele[0] as u8);
                found_alt += 1;
                //println!("allele {}" , allele[1] as char); // this will print the second base if the alt allele is for example AC (C).
            }
            if found_alt != 1 {
                panic!("Cannot manage multi-allelic SNPs!") // maybe simply skip?
            }
            let alte = get_sequence(&alt);
            // remember trim_alleles(&mut self)

            let rgenotypes = record.genotypes().unwrap();
            let genotypes = (0..self.reader.header.sample_count() as usize).map(|s| {
                                let geno_str = format!("{}", rgenotypes.get(s));
                                decode_genotype(geno_str)
                            }).collect_vec(); //useful to pre alloc?  Vec::<(bool, bool)>::with_capacity(self.samples.len());
 
            let pos = Coordinate { chr: chr, start: coord, end: coord+1}; // Right now only Snps.
            Some(Mutation { id: id, pos: pos, sequence_ref: refe, sequence_alt: alte, genotypes : genotypes})
        } else {
            return None;
        }
    }
}