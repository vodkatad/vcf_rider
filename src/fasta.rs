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
    pub sequence: Vec<u8>
}

pub struct FastaReader {
    reader: BufReader<File>,
    buffer: String
}

impl FastaReader {
    pub fn open(file: File) -> io::Result<FastaReader> {
        let mut b = String::new();
        let mut r = BufReader::new(file);
        return match r.read_line(&mut b) {
            Ok(_) => Ok(FastaReader { reader: r, buffer: b }),
            Err(e) => Err(e)
        }
    }

    pub fn open_path(path: &str) -> io::Result<FastaReader> {
        return match File::open(path) {
            Ok(file) => FastaReader::open(file),
            Err(e) => Err(e)
        };
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
        let mut sequence = String::new();
        let mut sequence_complete = false;
        while !sequence_complete && self.reader.read_line(&mut self.buffer).unwrap() > 0 {
            if self.buffer.chars().nth(0).unwrap() != '>' {
                sequence.push_str(self.buffer.trim_right());
                self.buffer.clear();
            }
            else {
                sequence_complete = true;
            }
        }
        return Some(Fasta { id: id, sequence: sequence.into_bytes() })
    }
}
