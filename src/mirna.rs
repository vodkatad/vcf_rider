// We need:
// Struct representing the miRNA (sequence, name and length)
// a miRNA reader
// implement CanScoreSequence for miRNAs.

// tab delimited file:
//data@tungsteno:/mnt/red/elly/bioinfotree/prj/roar/dataset/0.5/EUR/miRseeds_alteration$ head miRseeds_list 
//AAAAACC miR-548ac/548bb-3p/548d-3p/548h-3p/548z -1      hsa-miR-548ac
//AAAAACC miR-548ac/548bb-3p/548d-3p/548h-3p/548z -1      hsa-miR-548bb-3p
//AAAAACC miR-548ac/548bb-3p/548d-3p/548h-3p/548z -1      hsa-miR-548d-3p
//AAAAACC miR-548ac/548bb-3p/548d-3p/548h-3p/548z -1      hsa-miR-548h-3p
//AAAAACC miR-548ac/548bb-3p/548d-3p/548h-3p/548z -1      hsa-miR-548z
//AAAAACU miR-548aq-3p/548am-3p/548aj-3p/548ah-3p/548ae-3p/548j-3p/548x-3p        -1      hsa-miR-548ae-3p

// every subsequence score will be 1, 0:
// 1: has a seed of type 8mer, 7mer-A1 e 7mer-m8
// 0: doesn't have it
// 6mer are not considered as hits

// we always consider subsequence of length 8
// $motif = reverse $seed;
// $motif =~ tr/ACGT/TGCA/;
// $motif_6mer = substr $motif, 1, 6;
// $motif_m8 = substr $motif, 0, 1;
// $motif_A1 = 'A';

// $tentative_6mer_ref = substr $fasta_ref[$i], $j+1, 6;
// $tentative_m8_ref = substr $fasta_ref[$i], $j, 1;
// $tentative_A1_ref = substr $fasta_ref[$i], $j+7, 1;

/*        if ($tentative_6mer eq $motif_6mer)
        {
                $site = '6mer';
                $start = $j+1;
                $end = $j+6;
                if ($tentative_A1 eq $motif_A1)
                {
                        $site = '7mer-A1';
                        $end = $j+7;
                }
                if ($tentative_m8 eq $motif_m8 && $site eq '6mer')
                {
                        $site = '7mer-m8';
                        $start = $j;
                }
                if ($tentative_m8 eq $motif_m8 && $site eq '7mer-A1')
                {
                        $site = '8mer';
                        $start = $j;
                }
        }*/

/*
I get that I would need to check for == of substrings in this way:
- 7mer-A1: seq 1-6 (0based exclusive end) == to seed 6 last bases, seq 6==A
- 7mer-m8: seq 1-6 (0based exclusive end) == to seed 6 last bases, seq 0==first base of seed --> 0-7 ==seed
- 8mer: seq 1-6 (0based exclusive end) == to seed, seq 6==A, seq 0==first base of seed

So:
if seq == S1,S2,S3,S4,S5,S6,S7,A: 8mer
elseif seq start == S1,S2,S3,S4,S5,S6,S7: 7mer-m8
elseif seq end == S2,S3,S4,S5,S6,S7,A: 7mer-A1

better to compare whole string or first check the 6mer and then only single positions like the perl?

8mer if S1-S7A == my seq
7mer-A1 if S2-S7A == last 7 of my seq and first of my seq !=S1
7mer-m8 if S1-S7 == first 7 of my seed

Need to consider strand! Only the given one!
*/

use std::io;
use std::io::prelude::*;
use std::io::BufReader;
use std::fs::File;
use std::iter::Iterator;
use rider::CanScoreSequence;


const A: usize = 0;
const C: usize = 1;
const G: usize = 2;
const T: usize = 3;
const N: usize = 4;

/// Struct representing a miRNA seed
#[derive(Debug)]
pub struct Seed {
    /// The name of this miRNA family
    pub name: String,
    /// A `Vec<u8>` representing this miRNA seed
    pub sequence: Vec<u8>, // they will always have a length of 8, can we optimize? XXX TODO
}

// Implementing CanScoreSequence is the only step needed to use vcf_rider. 
// This is its implementation to count number of miRNA seed found on sequences 
// (we consider 8mer, 7mer-A1, and 7mer-m8 matches).
impl CanScoreSequence for Seed {
    fn get_length(&self) -> usize {
        self.sequence.len()
    }

    fn get_name(&self) -> &str {
        &self.name
    }

    /// Computes the score for this `Seed` on a given sequence starting at the given pos.
    /// The returned score is 1 if there is a match for this seed on the given sequence
    /// (8mer, 7mer-A1, and 7mer-m8 are considere match), 0 otherwise.
    /// The library only calls get score for sequences with lengths sufficient to host the seed starting at pos.
    /// # Arguments
    ///
    /// * `pos` - the position where we want to score this sequence
    /// * `sequence` - the score that we want to score
    fn get_score(&self, pos: usize, sequence: &[u8]) -> f64 {
        return 0f64;
        /*if seq == S1,S2,S3,S4,S5,S6,S7,A: 8mer
        elseif seq start == S1,S2,S3,S4,S5,S6,S7: 7mer-m8
        elseif seq end == S2,S3,S4,S5,S6,S7,A: 7mer-A1*/
    }
}

/// Struct used to read Seed. It will be implement an Iterator of `Seed` structs.
pub struct SeedReader {
    reader: BufReader<File>,
    buffer: String
}

impl SeedReader {
    /// Opens a tab delimited file with miRNA seeds file returning a `Result<SeedReader>`.
    /// # Arguments
    ///
    /// * `file` - the File with miRNA information
    ///
    /// # Errors 
    /// When the first line of the file cannot be read, the Error returned by read_line.
    /// Or the error returned by `File::open` if the file is not readable.
    pub fn open(file: File) -> io::Result<SeedReader> {
        let mut b = String::new();
        let mut r = BufReader::new(file);
        match r.read_line(&mut b) {
            Ok(_) => Ok(SeedReader { reader: r, buffer: b }),
            Err(e) => Err(e)
        }
    }

    /// Opens a tab delimited file with miRNA seeds file returning a `Result<SeedReader>`.
    /// # Arguments
    ///
    /// * `path` - the path to the miRNA file
    ///
    /// # Errors 
    /// When the first line of the file cannot be read, the Error returned by read_line.
    /// Or the error returned by `File::open` if the file is not readable.
    pub fn open_path(path: &str) -> io::Result<SeedReader> {
        match File::open(path) {
            Ok(file) => SeedReader::open(file),
            Err(e) => Err(e)
        }
    }
}

impl Iterator for SeedReader {
    type Item = Seed;

    fn next(&mut self) -> Option<Seed> {
        if self.buffer.is_empty() {
            return None;
        }
        let line = self.buffer.trim_right().to_owned();
        let mut tokens = line.split("\t");
        let name = tokens.nth(0).unwrap().to_owned();
        let mut sequence : Vec<u8> = Vec::with_capacity(8);
        for nuc in tokens.nth(0).unwrap().to_owned().as_bytes() {
            sequence.push(match *nuc {
                b'A' => 0u8,
                b'C' => 1u8,
                b'G' => 2u8,
                b'U' => 3u8,
                b'N' => 4u8,
                _ => panic!("Seed {} with a not allowed char {}", &name, *nuc as char),
            });
        }
        self.buffer.clear();
        self.reader.read_line(&mut self.buffer).unwrap();
        Some(Seed {name: name.to_owned(), sequence: sequence })
    }
}
