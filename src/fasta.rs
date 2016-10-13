use std::io;
use std::io::prelude::*;
use std::io::BufReader;
use std::fs::File;
use std::iter::Iterator;

/*struct fas {
    int *seq;
    char *id;
    double background[BASES];
    int length;
    int n_portions;
};
*/  
#[derive(Debug)]
pub struct Fasta {
    pub id: String,
    pub sequence: Vec<u8>,
    pub background: Vec<f64>
}

pub struct FastaReader {
    reader: BufReader<File>,
    buffer: String
}

impl FastaReader {
    pub fn open(file: File) -> io::Result<FastaReader> {
        let mut b = String::new();
        let mut r = BufReader::new(file);
        match r.read_line(&mut b) {
            Ok(_) => Ok(FastaReader { reader: r, buffer: b }),
            Err(e) => Err(e)
        }
    }

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
        let id = self.buffer.trim_right().to_owned();

        self.buffer.clear();

        // Read sequence lines and copy them to s.
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
        // TODO FIXME read bg from file/fasta
        Some(Fasta { id: id, sequence: sequence, background : vec!(0.298947240099661, 0.200854143743417, 0.200941012710477, 0.299257603446445)})
    }
}