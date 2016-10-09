use std::io;
use std::io::prelude::*;
use std::io::BufReader;
use std::fs::File;
use std::iter::Iterator;
use std::io::Lines;

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
    pub ll: Vec<Vec<f64>>,
    pub llrc: Vec<Vec<f64>>,
    pub freq: Vec<Vec<f64>>
}

pub struct MatrixReader {
    reader: BufReader<File>
}

impl MatrixReader {
    pub fn open(file: File) -> io::Result<MatrixReader> {
        let r = BufReader::new(file);
        return Ok(MatrixReader { reader: r });
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
        let mut name = "";
        let mut old_name = "";
        let mut freq : Vec<Vec<f64>> = Vec::new();
        let mut ll : Vec<Vec<f64>> = Vec::new();
        let mut llrc : Vec<Vec<f64>> =  Vec::new();;
        for line in self.reader.lines() {
            let mut tokens = line.unwrap().split("\t");
            name = tokens.nth(0).unwrap();
            if old_name != "" && name != old_name {
                // finished processing a matrix
                let result = Some(Matrix { name: old_name.to_owned(), ll: ll, llrc: llrc, freq: freq });
                freq = Vec::new();
                ll = Vec::new();
                llrc = Vec::new();
                return result;
            }
            let mut counts : Vec<f64> = Vec::with_capacity(BASES);
            for i in 0..BASES {
                let count: f64 = tokens.next().unwrap().parse::<f64>().unwrap();
                counts.push(count);
            }
            freq.push(counts);
            old_name = name;
        }
        return Some(Matrix {name: old_name.to_owned(), ll: ll, llrc: llrc, freq: freq });
    }
}
