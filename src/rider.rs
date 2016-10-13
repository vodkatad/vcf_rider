use bio::io::bed;
use std::fs;
use super::fasta;

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
    pub min_len: usize,
    pub max_len: usize,
    pub parameters: Vec<T>
    // TODO the operation to be used to manage scores
}

// TODO ->
pub fn get_scores<T : CanScoreSequence>(params: RiderParameters<T>, vcf_path: &str, mut bed_reader: bed::Reader<fs::File>, ref_path: &str) {
    println!("I would use {}", vcf_path);
    println!("With parameters {} {}", params.min_len, params.max_len);

    let referenceseq: fasta::Fasta = {
        if let Ok(mut reader) = fasta::FastaReader::open_path(ref_path) {
            let referenceseq = match reader.next() {
                Some(x) => x,
                None => panic!("Empty fasta?")
            };
            let other = reader.next(); // Very inefficient to read all of it, we will read two fasta before giving the error.
            match other {
                Some(x) => x,
                None => panic!("Right now you can use this only on single chr! You have given either a multifasta")
            };
            referenceseq
        } else {
            panic!("Error reading reference fasta {}")
        }
    };

    let used_chr = referenceseq.id;
    println!("Fasta ref {}", used_chr);
    // load vcf -> open file, skip # headers, first real entry
    // We could use a VcfReader similar to others.

    // chr check 
    for r in bed_reader.records() {
        let record = r.ok().expect("Error reading record");
        println!("bed name: {}", record.name().expect("Error reading name"));
        // chr check

        // let mut pos = r.start
        // while pos + params.max_len < r.end
            // obtain_seq(r.chr, r.start, params.max_len, VcfReader)
            // this will give us 2^n seqs where n in the n of snps found in r.start-r.rstart+params.max_len 
            // seqs will be ordered in a specific order: the first one is the reference one and TODO
            // for every p params.parameters
                // p.get_score(0usize, seq)
                // this will give us a vector  

        // find_overlapping_snps(r.chr, r.start, r.end, VcfReader, snps_on_bed)
        // We assume to receive ordered vcf and bed therefore we can skip vcf entries < bed,
        // if vcf > bed break and go to next bed
        // otherwise we have to renew the given snps_on_bed parameter according to bed coords if we decide to allow bed overlaps.
        // Without overlapping bed snp_on_bed param is not needed, it will always be renewed

    }   
}