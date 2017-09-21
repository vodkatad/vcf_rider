use std::io;
use std::io::prelude::*;
use std::io::BufReader;
use std::fs::File;
use std::iter::Iterator;

/// Struct representing a fasta: its id, then the sequence (encoded by a vector of u8, with ACGTN -> 01234) and the background frequencies of nucleotides (needed
/// to compute TBA values)
#[derive(Debug)]
pub struct Fasta {
    /// The id of this fasta
    pub id: String,
    /// The u8 encoded sequence
    pub sequence: Vec<u8>,
    /// Background frequences of ACGT in this fasta
    pub background: Vec<f64>
}

/// Struct used to read fasta files. It will be implement an Iterator of Fasta structs.
/// The String buffer is used to get characters from the BufReader one by one and convert
/// them to our internal Vec<u8> representation.
///
/// # Panics 
/// While iterating on it, if there is a not allowed nucleotide (only ACTGN are allowed).
pub struct FastaReader {
    reader: BufReader<File>,
    buffer: String
}

impl FastaReader {
    /// Opens a fasta file returning a Result<FastaReader>.
    /// # Arguments
    ///
    /// * `file` - the fasta File
    ///
    /// # Errors 
    /// When the first file of the fasta file cannot be read, the Error returned by read_line.
    pub fn open(file: File) -> io::Result<FastaReader> {
        let mut b = String::new();
        let mut r = BufReader::new(file);
        match r.read_line(&mut b) {
            Ok(_) => Ok(FastaReader { reader: r, buffer: b }),
            Err(e) => Err(e)
        }
    }
    
    /// Opens a fasta file returning a Result<FastaReader>.
    /// # Arguments
    ///
    /// * `path` - the path to the fasta file
    ///
    /// # Errors 
    /// When the first file of the fasta file cannot be read, the Error returned by read_line
    pub fn open_path(path: &str) -> io::Result<FastaReader> {
        match File::open(path) {
            Ok(file) => FastaReader::open(file),
            Err(e) => Err(e)
        }
    }
}

impl Iterator for FastaReader {
    type Item = Fasta;

    fn next(&mut self) -> Option<Fasta> {
        if self.buffer.is_empty() {
            return None;
        }
        // The first line is loaded in self.buffer when the FastaReader is created (and when finishing a fasta for multifasta files).
        let id = self.buffer.trim_right().to_owned();

        self.buffer.clear();

        // Read sequence lines and copy them to sequence with our u8 encoding for nucleotides.
        let mut sequence : Vec<u8> = Vec::new();
        while self.reader.read_line(&mut self.buffer).unwrap() > 0 {
            if self.buffer.chars().nth(0).unwrap() != '>' {
                for nuc in self.buffer.trim_right().as_bytes() {
                    sequence.push(match *nuc {
                        b'A' => 0u8,
                        b'C' => 1u8,
                        b'G' => 2u8,
                        b'T' => 3u8,
                        b'N' => 4u8,
                        _ => panic!("Fasta {} with a not allowed char {}", &id, *nuc as char),
                    });
                }
                self.buffer.clear();
            }
            else {
                break;
            }
        }
        // TODO FIXME compute bg. Right now we have hard encoded the one used in all our project (historic hg19 intergenic by Ivan).
        Some(Fasta { id: id, sequence: sequence, background : vec!(0.298947240099661, 0.200854143743417, 0.200941012710477, 0.299257603446445)})
    }
}