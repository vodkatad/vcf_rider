extern crate bio;

use bio::io::bed;
use std::iter::MinMaxResult::{self, NoElements, OneElement, MinMax};

/// Our vcf_rider main function will receive a Vec<T: CanScoreSequence>
/// and call it for every T on subsequences of the genomes of the samples 
/// doing it only for each variant subsequence once. 
/// This trait will need to be able to compute a score on a given sequence,
/// represented by a splice of an array of u8 [TODO] starting for a given
/// position (it is guaranteed by the lib that the used position will be given
/// inside the sequence, i.e. sequence.len() - self.get_length() >= 0).
pub trait CanScoreSequence {
    /// Returns a score for the given sequence starting at position pos.
    ///
    /// # Arguments
    ///
    /// * `self` - the object with trait CanScoreSequence.
    /// * `pos` - the position in the sequence where the score will be calculated.
    ///           The check that sequence.len() - self.get_length() >= 0 IS NOT DONE HERE. 
    /// * `sequence`- the sequence that needs to be scored, encoded as [ACGTN]-[01234]
    fn get_score(&self, pos: usize, sequence: &[u8]) -> f64;

    /// Returns the length of sequence that this object can score.
    ///
    /// # Arguments
    ///
    /// * `self` - the object with trait CanScoreSequence.
    fn get_length(&self) -> usize;
}

pub struct RiderParameters<T : CanScoreSequence> {
    min_max_len: MinMaxResult,
    parameters: Vec<T>
    // TODO the operation to be used to manage scores

}

// TODO -> (but without an iterator it would not be lazy?)
pub fn get_scores<T : CanScoreSequence>(params: RiderParameters<T>, vcf_path: &str, bed::Reader: bed_reader) {
    println!("I would use {}", vcf_path);
    println!("With parameters {:?}", params.min_max_len);
    for r in bed_reader.records() {
        let record = r.ok().expect("Error reading record");
        println!("bed name: {}", record.name().expect("Error reading name"));
    }   
}