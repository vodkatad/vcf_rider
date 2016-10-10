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

// The matrix is kept rows by columns, i.e. elements at indices 0-3 are the
// first row, 4-7 the second and so on. This makes it easier to append rows
// while reading them from a file.

#[derive(Debug)]
pub struct Matrix {
    rows: Vec<f64>
}

impl Matrix {
    pub fn new() -> Matrix {
        Matrix { rows: Vec::<f64>::new() }
    }
    
    pub fn with_capacity(len: usize) -> Matrix {
        Matrix { rows: Vec::<f64>::with_capacity(len) }
    }

    pub fn len(&self) -> usize {
        (self.rows.len() / BASES) as usize
    }

    pub fn push_row(&mut self, a: f64, c: f64, t: f64, g: f64) {
        self.rows.push(a);
        self.rows.push(c);
        self.rows.push(t);
        self.rows.push(g);
    }

    #[inline]
    pub fn get(&self, row: usize, col: usize) -> f64 {
        self.rows[row * BASES + col]
    }
}


#[derive(Debug)]
pub struct PWM {
    pub name: String,
    pub ll: Matrix,
    pub llrc: Matrix,
    pub freq: Matrix
}

pub struct PWMReader {
    reader: BufReader<File>,
    buffer: String
}

impl PWMReader {
    pub fn open(file: File) -> io::Result<PWMReader> {
        let mut b = String::new();
        let mut r = BufReader::new(file);
        match r.read_line(&mut b) {
            Ok(_) => Ok(PWMReader { reader: r, buffer: b }),
            Err(e) => Err(e)
        }
    }

    pub fn open_path(path: &str) -> io::Result<PWMReader> {
        match File::open(path) {
            Ok(file) => PWMReader::open(file),
            Err(e) => Err(e)
        }
    }
}

impl Iterator for PWMReader {
    type Item = PWM;

    fn next(&mut self) -> Option<PWM> {
        if self.buffer.is_empty() {
            return None;
        }
        let mut name;
        let mut old_name = "".to_owned();
        let mut freq = Matrix::new();
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
                let a: f64 = tokens.next().unwrap().parse::<f64>().unwrap();     
                let c: f64 = tokens.next().unwrap().parse::<f64>().unwrap();
                let t: f64 = tokens.next().unwrap().parse::<f64>().unwrap();
                let g: f64 = tokens.next().unwrap().parse::<f64>().unwrap();
                freq.push_row(a, c, t, g);
                old_name = name;
            }
            if self.reader.read_line(&mut self.buffer).unwrap() == 0 {
                finished = true;
                self.buffer.clear();
            }
        }
        let ll = Matrix::with_capacity(freq.len());
        let llrc =  Matrix::with_capacity(freq.len());;
        Some(PWM {name: old_name.to_owned(), ll: ll, llrc: llrc, freq: freq })
    }
}
