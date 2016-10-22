use bio::io::bed;
use std::fs;
use super::fasta;
use super::mutations;

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

    /// Returns the name of this.
    ///
    /// # Arguments
    ///
    /// * `self` - the object with trait CanScoreSequence.
    fn get_name(&self) -> &str;
}

pub struct RiderParameters<'a, T : CanScoreSequence + 'a> {
    pub min_len: usize,
    pub max_len: usize,
    pub parameters: &'a Vec<T>
    // TODO the operation to be used to manage scores http://doc.rust-lang.org/core/ops/
}

// TODO ->
pub fn get_scores<T : CanScoreSequence>(params: RiderParameters<T>, vcf_path: &str, mut bed_reader: bed::Reader<fs::File>, ref_path: &str) {

    // #[cfg(debug_assertions)]   attributes on non-item statements and expressions are experimental. (see issue #15701)
    //{
    println!("I would use {}", vcf_path);
    println!("With parameters {} {}", params.min_len, params.max_len);
    //}
    let referenceseq: fasta::Fasta = {
        if let Ok(mut reader) = fasta::FastaReader::open_path(ref_path) {
            let referenceseq = match reader.next() {
                Some(x) => x,
                None => panic!("Empty fasta?")
            };
            let other = reader.next(); // Very inefficient to read all of it, we will read two fasta before giving the error.
            match other {
                None => (),
                Some(_) => panic!("Right now you can use this only on single chr! You have given either a multifasta")
            };
            referenceseq
        } else {
            panic!("Error reading reference fasta")
        }
    };

    let used_chr = referenceseq.id;
    println!("Fasta ref {}", used_chr);
    // load vcf -> open file, skip # headers, first real entry
    // We could use a VcfReader similar to others.
    if let Ok(vcf) = mutations::VcfReader::open_path(vcf_path) {
        for sample in & vcf.samples {
            println!("sample {}", sample);
        }
        for snp in vcf {
            println!("snp {} {:?} {}", snp.coord, snp.sequence_ref, snp.id);
        }
    }

    // chr check
    for r in bed_reader.records() {
        let record = r.ok().expect("Error reading record");
        println!("bed name: {}", record.name().expect("Error reading name"));
        // chr check

        // let mut pos = r.start
        // while pos < r.end   // we do not do pos + params.max_len < r.end to avoid cumbersome management for the last portion
            // obtain_seq(r.chr, r.start, params.max_len, VcfReader, snps_on_seq)
            // this will give us 2^n seqs where n in the n of snps found in r.start-r.rstart+params.max_len
            // seqs will be ordered in a specific order: the first one is the reference one and the last one
            // the one with all mutated alleles.  Every SNP status is encoded by 0 if it is reference and 1 if it is
            // mutated. The first snp in the seq is encoded by the least significant position.
            // In this way it will be easy to build for each individual chr the indexes linking
            // them to their sequences. There will be a vector of tuples (usize,usize) for this.
            // obtain seq will return sequences of length params.max_len if possible otherwise shorter
            // ones and we will check to call get_score only on the right parameters
            // for every p params.parameters call on seq
            for i in 0..(*params.parameters).len() {
                println!("pwm {}", (*params.parameters).get(i).unwrap().get_name());
            }
                // if  p.len() < seq.len // not cumbersome but inefficient?
                // for every index present in the two vectors:
                    // p.get_score(0usize, seq)
                    // now we need to sum (or smt else) the scores assigning them to the right individuals.


    }
}

// find_overlapping_snps(r.start, r.end, VcfReader, snps_on_seq) - to pass coords use bed struct or a struct by me to easily add chr.
// We assume to receive ordered vcf and bed therefore we can skip vcf entries < r.start,
// if vcf > r.start+r.end empty snps_on_seq and return ---> not empty! leave them there but return false to avoid loosing vcf info
// otherwise we have to renew the given snps_on_seq parameter according to coords
// return true if there are overlapping snps and false otherwise

// obtain_seq(r.chr, r.start, params.max_len, VcfReader, snps_on_seq)
// snps_on_seq will be empty or contain snps found in the previous window
// we will call find_overlapping_snps
// if it returns false we get the reference sequence and return only it --- we need a struct for the return type of obtain_seq
// otherwise we need to build the sequences and the individual vectors

// Needed structs:
// coords and comparisons methods
// return type of obtain_seq with sequences and individual info
// return type of get scores with bed info, n of snps, [len of individuals seqs], scores for both alleles for individuals and individual ids
// vcf entry ?