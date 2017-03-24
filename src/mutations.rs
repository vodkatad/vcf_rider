use rust_htslib::bcf;
use std::io;
use itertools::Itertools;
use std::path::Path;

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
    pub pos: Coordinate, // Always end=start+1 to compute overlaps safely.
    pub sequence_ref: Vec<u8>,
    pub sequence_alt: Vec<u8>,
    pub genotypes: Vec<(bool, bool)>,
    pub is_indel: bool,
    pub indel_len: i64 // We use this to store the length of indels instead of pos.end.
    // the length of insertions is negative
}

pub struct VcfReader {
    reader: bcf::Reader,
    accept_phased: bool,
    pub samples: Vec<String>
}

// Adding a layer of abstration, I am not sure that we will use the lib.
impl VcfReader {
    pub fn open_path(path: &str, accept_phased: bool) -> io::Result<VcfReader> {
        match bcf::Reader::from_path(Path::new(path)) {
            Ok(reader) => {
                let samples = reader.header.samples().into_iter().map(|sample| String::from_utf8(sample.to_owned()).unwrap()).collect();
                Ok(VcfReader { reader: reader, accept_phased: accept_phased, samples: samples })
            }
            Err(_) => Err(io::Error::last_os_error())
            // How do errors work? rust_htslib::bcf::BCFError
        }
    }
}

fn get_sequence(seq : &Vec<u8>) -> Vec<u8> {
    let mut res : Vec<u8> = Vec::<u8>::with_capacity(1);
    if seq.iter().map(|x| *x as char).join("") == "<DEL>" {
        return vec![6u8, 6u8, 6u8]; 
        // Funny me. Marking big deletions as the number of the beast. Their length is needed to mark them as indels.
        // iter here and then to match nucleotides to numbers is bad?
    }
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
    /*if res.len() > 1 {
        panic!("Right now we do not handle more than SNPs!");
    }*/
    res
}
fn decode_allele(c: char) -> bool {
    match c {
        '0' => false,
        '1' => true,
         x => panic!("Wrongly encoded snp in vcf {:?}", x)
    }
}

fn decode_genotype(geno: String, accept_phased: bool) -> (bool, bool) {
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
    if sep != '|' && !accept_phased {
        panic!("Only phased genotypes are supported! {}", sep)
    }
    //println!("allele1 {:?}", a1);
    //println!("allele1 {:?}", a2);
    (decode_allele(a1),decode_allele(a2))
}

impl Iterator for VcfReader {
    type Item = Mutation;

    fn next(&mut self) -> Option<Mutation> {
        let mut record : bcf::Record = bcf::Record::new();
        if let Ok(_) = self.reader.read(&mut record) {
            let id = "".to_owned();
            let chr = "?".to_owned(); // where is it? rid?
            let mut coord = record.pos() as u64; // Manually checked: the lib converts 1 based coords of vcf to 0 based. bed records have u64
            let alleles = record.alleles().into_iter().map(|a| a.to_owned()).collect_vec(); // I do not like this to_owned...
            let mut alt  : Vec<u8> = Vec::<u8>::with_capacity(1);
            let mut found_alt = 0;
            let refe : Vec<u8> = get_sequence(&alleles[0]);
            //println!("ref {}", alleles[0].iter().map(|x| *x as char).join(""));
            for allele in alleles[1..].iter() { // why skip the first? The first is the reference. Then are listed all the alternative ones.
                alt.extend(allele.iter().map(|x| *x as u8)); // not the best way, play with scopes and mut and so on TODO
                found_alt += 1;
                //println!("allele {}" , allele[1] as char); // this will print the second base if the alt allele is for example AC (C).
            }
            if found_alt != 1 {
                panic!("Cannot manage multi-allelic SNPs! {:?}", alt) // maybe simply skip?
            }
            let mut alte = get_sequence(&alt);
            let mut len = 1;
            let indel = refe.len() != 1 || alte.len() != 1;
            if alte == vec![6u8, 6u8, 6u8] {
                let info_end = record.info("END".as_bytes()).integer();
                //let end = match info_end.unwrap_or_else(panic!("VCF with a <DEL> and no END info!")) {let end = match info_end.unwrap_or_else(panic!("VCF with a <DEL> and no END info!")) {
                let end = match info_end.unwrap() {
                    Some(l) => l,
                    None => panic!("VCF with a <DEL> and an END info not good"),
                };
                if end.len() > 1 || end[0] < 0 {
                    panic!("Multiallelic snp or wrong end");
                }
                len = (end[0] as i64 -1) - coord as i64; // minus one to go back to 0 based coords even for this one (I believe that the lib does not fix this).
                // Following vcf v4.2 specs the length is only "approximate", but it is ok to use end - pos (end is exclusive).
            } else if indel {
                let ref_len = refe.len();
                let alt_len = alte.len(); // what is the cost of calling .len()?
                len = (ref_len as i64 - alt_len as i64) as i64;
                // We add 1 to the starting coord because indel are encoded with the first base 
                // of reference out of the insertion or deletion.
                coord += 1;
                // and we remove the first base from the alt allele
                alte.remove(0);
                // IN have a negative length, DEL a positive one.
            }
            let pos = Coordinate { chr: chr, start: coord, end: coord+1}; 
            //println!("alt {}", alt.iter().map(|x| *x as char).join(""));
            // remember trim_alleles(&mut self)

            let rgenotypes = record.genotypes().unwrap();
            let genotypes = (0..self.reader.header.sample_count() as usize).map(|s| {
                                let geno_str = format!("{}", rgenotypes.get(s));
                                decode_genotype(geno_str, self.accept_phased)
                            }).collect_vec(); //useful to pre alloc?  Vec::<(bool, bool)>::with_capacity(self.samples.len());
 
            println!("ref {:?} alt {:?} s {} e{} len {}", refe, alte, pos.start, pos.end, len);
            Some(Mutation { id: id, pos: pos, sequence_ref: refe, sequence_alt: alte, genotypes : genotypes, is_indel : indel, indel_len : len})
        } else {
            return None;
        }
    }
}
