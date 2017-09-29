use std::io;
use std::io::prelude::*;
use std::io::BufReader;
use std::fs::File;
use std::iter::Iterator;
use rider::CanScoreSequence;

const BASES: usize = 4;
const NN: f64 = 0f64;
const A: usize = 0;
const C: usize = 1;
const G: usize = 2;
const T: usize = 3;
const N: usize = 4;

/// Struct used to represent PWM likelihoods/nucleotide frequencies. 
/// ncols is used to be able to represent 4 and 5 nucleotides PWMs (to manage N in sequences).
/// No consistency checks are performed for get/set operations for efficiency.
#[derive(Debug)]
pub struct Matrix {
    // The matrix is kept rows by columns, i.e. elements at indices 0-3 are the
    // first row, 4-7 the second and so on. This makes it easier to append rows
    // while reading them from a file.
    /// Vector that contains the linearized values of the matrix, row by row.
    rows: Vec<f64>,
    /// Number of columns of this `Matrix`
    ncols: usize
}

impl Matrix {
    /// Creates a new matrix with n columns.
    /// # Arguments
    ///
    /// * `n` - the number of columns
    pub fn new(n: usize) -> Matrix {
        Matrix { rows: Vec::<f64>::new(), ncols: n }
    }
    
    /// Creates a new matrix preallocated for len elements.
    /// Warning: len is not checked to be compatible with n.
    ///
    /// # Arguments
    ///
    /// * `n` - the number of columns
    /// * `len` - the number of elements to be preallocated
    pub fn with_capacity(n: usize, len: usize) -> Matrix {
        Matrix { rows: Vec::<f64>::with_capacity(len), ncols: n }
        // Here instead a check on len % n could be done, maybe.
    }

    #[inline]
    /// Returns the number of rows of this `Matrix`
    pub fn len(&self) -> usize {
        // TODO explain why inline. https://internals.rust-lang.org/t/when-should-i-use-inline/598
        (self.rows.len() / self.ncols) as usize
    }

    /// Adds a 4 element row to this `Matrix`. No checks of compatibility with `self.ncols`.
    ///
    /// # Arguments
    ///
    /// * `a` - the likelihood for nucleotide A that we want to add
    /// * `c` - the likelihood for nucleotide C that we want to add
    /// * `g` - the likelihood for nucleotide G that we want to add
    /// * `t` - the likelihood for nucleotide T that we want to add
    pub fn push_row_4(&mut self, a: f64, c: f64, g: f64, t: f64) {
        // We do not do an if to check if we are breaking up
        // our 2d matrix calling the wrong function for efficiency's sake.
        // TODO FIXME generalize.
        self.rows.push(a);
        self.rows.push(c);
        self.rows.push(g);
        self.rows.push(t);
    }

    /// Adds a 5 element row to this `Matrix`. No checks of compatibility with `self.ncols`.
    ///
    /// # Arguments
    ///
    /// * `a` - the likelihood for nucleotide A that we want to add
    /// * `c` - the likelihood for nucleotide C that we want to add
    /// * `g` - the likelihood for nucleotide G that we want to add
    /// * `t` - the likelihood for nucleotide T that we want to add
    /// * `n` - the likelihood for nucleotide N that we want to add
    pub fn push_row_5(&mut self, a: f64, c: f64, g: f64, t: f64, n: f64) {
        self.rows.push(a);
        self.rows.push(c);
        self.rows.push(g);
        self.rows.push(t);
        self.rows.push(n);
    }

    #[inline]
    /// Returns the element in row, col of this `Matrix`. No boundary checks are performed!
    pub fn get(&self, row: usize, col: usize) -> f64 {
        self.rows[row * self.ncols + col]
    }
}

/// Struct representing a whole Positional Weight Matrix.
#[derive(Debug)]
pub struct PWM {
    /// The name of this PWM
    pub name: String, // &str here?
    /// A `Matrix` storing likelihoods for this PWM
    pub ll: Matrix,
    /// A `Matrix` storing reverse complement likelihoods for this PWM. i.e. its elements 0,0 has the likelihood
    /// of a T (reverse of A) in the last position of the PWM - used to efficiently score the reverse complement sequence while
    /// analyzing sequences going forward.
    pub llrc: Matrix,
    /// A `Matrix` storing the original frequencies for nucleotides, these becomes likelihoods when divided by a given sequence background.
    pub freq: Matrix
}

// Implementing CanScoreSequence is the only step needed to use vcf_rider. 
// This is its implementation for PWM to compute Total Binding Affinities.
impl CanScoreSequence for PWM {
    fn get_length(&self) -> usize {
        self.freq.len()
    }

    fn get_name(&self) -> &str {
        &self.name
    }

    // pub here is unnecessary? Why? FIXME all pub/no pub needs fixing/review.

    /// Computes the score for this `PWM` on a given sequence starting at the given pos.
    /// The returned score is the highest between the "straight" one and the one on the reverse complement.
    /// We do not need checks on lengths here because we will
    /// perform then in the library at a (hopefully) least repeated step.
    /// The library only calls get score for sequences with lengths sufficient to host the PWM starting at pos.
    /// # Arguments
    ///
    /// * `pos` - the position where we want to score this sequence
    /// * `sequence` - the score that we want to score
    fn get_score(&self, pos: usize, sequence: &[u8]) -> f64 {
        let mut res1: f64 = 1f64;
        let mut res2: f64 = 1f64;
        let mut j = pos;
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

// Implementing all the functions needed to load and then use a PWM.
impl PWM {
    /// Given a `PWM` that has frequencies loaded and some background frequencies it fills `ll` and ``ll_rc``.
    /// # Arguments
    ///
    /// * `bg` - a vector with A C G and T nucleotide frequencies to be considered as background for this PWM.
    pub fn compute_ll(&mut self, bg: &Vec<f64>) {
        let mut frac : Vec<f64> = Vec::with_capacity(4); // TODO fire up tests and change to 5. BASES+1
        let matrix_len = self.freq.len(); // useful? Should not be due to automatic inlining.
    	for j in 0..matrix_len {
	    	for i in 0..BASES {
                frac.push(self.freq.get(j, i) / bg[i]);
                // Corresponding C (log2_ratio there is performed for some combo of the binary options):
                //m->ll[j][i] = log2_ratio(m->freq[j][i], bg[i], info);

		    }
            // We always score N as 0. Every sequence that we will score where an N will appear in the scored
            // subsequence will have 0 since we are multipling single positions scores.
		    frac.push(NN);
            self.ll.push_row_5(frac[0], frac[1], frac[2], frac[3], frac[4]);
            frac.clear();
	    }
    	for j in 0..matrix_len {
	    	for i in 0..BASES {
                // We are now getting reverse complement positions from the PWM, so starting from the end for rows (matrix_len-0-1= last row)
                // and accessing complement columns. We use the ll that we have filled in the previous loop.
                frac.push(self.ll.get(matrix_len-j-1, match i {
                    A => T,
                    C => G,
                    G => C,
                    T => A,
                    N => N,
                    _ => panic!("The programmer this time screwed up big time!")
                }))
                // Corresponding C:
			    //m->llrc[j][i] = m->ll[(m->length) - j - 1][encoded_rc(i)];
		    }
		    frac.push(NN);
            //m->llrc[j][N] = NN;
            self.llrc.push_row_5(frac[0], frac[1], frac[2], frac[3], frac[4]);
            frac.clear();
	    }
    }
}

/// Struct used to read PWM. It will be implement an Iterator of `PWM` structs.
/// The String buffer is used to get characters from the `BufReader` one by one and load
/// the `freq` for the loaded `PWM` structs.
/// The expected format is the same used by matrix_rider in C, therefore a tab delimited file
/// with these columns:
/// pwm_name, pos, count_a, count_c, count_g, count_t
/// pos is ignored, count_x should be positive integer different from 0.
pub struct PWMReader {
    reader: BufReader<File>,
    buffer: String
}

impl PWMReader {
    /// Opens a PWM file returning a `Result<PWMReader>`.
    /// # Arguments
    ///
    /// * `file` - the PWM File
    ///
    /// # Errors 
    /// When the first file of the file cannot be read, the Error returned by read_line.
    pub fn open(file: File) -> io::Result<PWMReader> {
        let mut b = String::new();
        let mut r = BufReader::new(file);
        match r.read_line(&mut b) {
            Ok(_) => Ok(PWMReader { reader: r, buffer: b }),
            Err(e) => Err(e)
        }
    }

    /// Opens a pwm file returning a `Result<PWMReader>`.
    /// # Arguments
    ///
    /// * `path` - the path to the PWM file
    ///
    /// # Errors 
    /// When the first file of the file cannot be read, the Error returned by read_line
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
        let mut freq = Matrix::new(BASES);
        let mut finished = false;
        // We have finished a PWM when the name changes (first element of our rows).
        while !finished {
            let line = self.buffer.trim_right().to_owned();
            let mut tokens = line.split("\t");
            name = tokens.nth(0).unwrap().to_owned();
            tokens.next(); // We skip the second column with position.
            if old_name != "" && name != old_name {
                // finished processing a matrix
                finished = true;
            } else {
                let a: u64 = tokens.next().unwrap().parse::<u64>().unwrap();     
                let c: u64 = tokens.next().unwrap().parse::<u64>().unwrap();
                let g: u64 = tokens.next().unwrap().parse::<u64>().unwrap();
                let t: u64 = tokens.next().unwrap().parse::<u64>().unwrap();
                let tot: f64 = (a + c + g + t) as f64;
                // TODO ADD error checking, this is get_fraction_from_pcounts FIXME.
                // No zeroes etc etc.
                freq.push_row_4(a as f64/ tot, c  as f64/ tot, g  as f64/ tot, t  as f64/ tot);
                old_name = name;
            }
            if !finished {
                self.buffer.clear();
                if self.reader.read_line(&mut self.buffer).unwrap() == 0 {
                    finished = true;
                    self.buffer.clear();
                }
            }
        }
        let ll = Matrix::with_capacity(BASES+1, freq.rows.len());
        let llrc =  Matrix::with_capacity(BASES+1, freq.rows.len());;
        Some(PWM {name: old_name.to_owned(), ll: ll, llrc: llrc, freq: freq })
    }
}
