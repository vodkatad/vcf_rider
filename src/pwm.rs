use std::io;
use std::io::prelude::*;
use std::io::BufReader;
use std::fs::File;
use std::iter::Iterator;

/*struct matrix_ll_ {
    double **ll;
    double **llrc;
    double **freq;
    int length;
    char *name;
};*/
const BASES: usize = 4;

#[derive(Debug)]
pub struct Matrix {
    pub name: String,
    pub ll: Vec<f64>,
    pub llrc: Vec<f64>,
    pub freq: Vec<f64>
}

pub struct MatrixReader {
    reader: BufReader<File>,
    buffer: String
}

impl MatrixReader {
    pub fn open(file: File) -> io::Result<MatrixReader> {
        let mut b = String::new();
        let mut r = BufReader::new(file);
        return match r.read_line(&mut b) {
            Ok(_) => Ok(MatrixReader { reader: r, buffer: b }),
            Err(e) => Err(e)
        }
    }

    pub fn open_path(path: &str) -> io::Result<MatrixReader> {
        return match File::open(path) {
            Ok(file) => MatrixReader::open(file),
            Err(e) => Err(e)
        };
    }
}

impl Iterator for MatrixReader {
    type Item = Matrix;

    fn next(&mut self) -> Option<Matrix> {
        if self.buffer.is_empty() {
            return None;
        }
        let mut name;
        let mut old_name = "".to_owned();
        let mut freq : Vec<f64> = Vec::new();
        let mut finished = false;
        while !finished {
            let line = self.buffer.trim_right().to_owned();
            self.buffer.clear();
            let mut tokens = line.split("\t");
            name = tokens.nth(0).unwrap().to_owned();
            tokens.next(); // We skip the second column with position.
            if old_name != "" && name != old_name {
                // finished processing a matrix
                finished = true;
            } else {
                for _i in 0..BASES {
                    let count: f64 = tokens.next().unwrap().parse::<f64>().unwrap();
                    freq.push(count);
                }
                old_name = name;
            }
            if self.reader.read_line(&mut self.buffer).unwrap() == 0 {
                finished = true;
                self.buffer.clear();
            }
        }
        let ll : Vec<f64> = Vec::with_capacity(freq.len());
        let llrc : Vec<f64> =  Vec::with_capacity(freq.len());;
        Some(Matrix {name: old_name.to_owned(), ll: ll, llrc: llrc, freq: freq })
    }
}
