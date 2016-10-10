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
const NN: f64 = 0f64;
const A: usize = 0;
const C: usize = 1;
const G: usize = 2;
const T: usize = 3;
const N: usize = 4;

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

    #[inline]
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

//TODO Matrix<usize> -> usize is NCOL.

#[derive(Debug)]
pub struct MatrixN {
    rows: Vec<f64>
}

impl MatrixN {
    pub fn new() -> MatrixN {
        MatrixN { rows: Vec::<f64>::new() }
    }
    
    pub fn with_capacity(len: usize) -> MatrixN {
        MatrixN { rows: Vec::<f64>::with_capacity(len) }
    }

    #[inline]
    pub fn len(&self) -> usize {
        (self.rows.len() / (BASES+1)) as usize
    }

    pub fn push_row(&mut self, a: f64, c: f64, t: f64, g: f64, n: f64) {
        self.rows.push(a);
        self.rows.push(c);
        self.rows.push(t);
        self.rows.push(g);
        self.rows.push(n);
    }

    #[inline]
    pub fn get(&self, row: usize, col: usize) -> f64 {
        self.rows[row * (BASES+1) + col]
    }
}



#[derive(Debug)]
pub struct PWM {
    pub name: String,
    pub ll: MatrixN,
    pub llrc: MatrixN,
    pub freq: Matrix
}

impl PWM {
    pub fn compute_ll(&mut self, bg: &Vec<f64>) {
        let mut frac : Vec<f64> = Vec::with_capacity(4);
        let matrix_len = self.freq.len(); // useful?
    	for j in 0..matrix_len {
	    	for i in 0..BASES {
                frac.push(self.freq.get(j, i) / bg[i]);
                //m->ll[j][i] = log2_ratio(m->freq[j][i], bg[i], info);

		    }
		    frac.push(NN);
            self.ll.push_row(frac[0], frac[1], frac[2], frac[3], frac[4]);
            frac.clear();
	    }
    	for j in 0..matrix_len {
	    	for i in 0..BASES {
                frac.push(self.ll.get(matrix_len-j-1, match i {
                    A => T,
                    C => G,
                    G => C,
                    T => A,
                    N => N,
                    _ => panic!("The programmer this time screwed up big time!")
                }))
			    //m->llrc[j][i] = m->ll[(m->length) - j - 1][encoded_rc(i)];
		    }
		    frac.push(NN);
            //m->llrc[j][N] = NN;
            self.llrc.push_row(frac[0], frac[1], frac[2], frac[3], frac[4]);
            frac.clear();
	    }
    }

    pub fn get_affinity(&self, sequence: &[u8]) -> f64 {
        let mut res1: f64 = 1f64;
        let mut res2: f64 = 1f64;
        let mut j = 0;
        for i in 0..self.freq.len() {
            res1 *= self.ll.get(i, sequence[j] as usize);
            res2 *= self.llrc.get(i, sequence[j] as usize);
            j += 1;
        }
        if res1 > res2 {
            res1
        }
        else {
            res2
        }
    }
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
                let a: u64 = tokens.next().unwrap().parse::<u64>().unwrap();     
                let c: u64 = tokens.next().unwrap().parse::<u64>().unwrap();
                let t: u64 = tokens.next().unwrap().parse::<u64>().unwrap();
                let g: u64 = tokens.next().unwrap().parse::<u64>().unwrap();
                let tot: f64 = (a + c + t + g) as f64;
                // TODO ADD error checking, this is get_fraction_from_pcounts
                freq.push_row(a as f64/ tot, c  as f64/ tot, t  as f64/ tot, g  as f64/ tot);
                old_name = name;
            }
            if self.reader.read_line(&mut self.buffer).unwrap() == 0 {
                finished = true;
                self.buffer.clear();
            }
        }
        let ll = MatrixN::with_capacity(freq.rows.len());
        let llrc =  MatrixN::with_capacity(freq.rows.len());;
        Some(PWM {name: old_name.to_owned(), ll: ll, llrc: llrc, freq: freq })
    }
}
