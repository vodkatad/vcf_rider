use bio::io::bed;
use std::fs;
use super::fasta;
use super::mutations;
use std::collections::VecDeque;

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
            println!("snp {:?} {:?} {}", snp.pos, snp.sequence_ref, snp.id);
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

/// Function that advances on the VcfReader (Iterator of Mutation) until the first snp that does not overlap with the given window, putting
/// in snps_buffer all the overlapping snps and their number and then the first not overlapping snp.
///
/// # Arguments
///
/// * `window` - the window where we search for overlapping SNPs. This function should be called giving them in order.
/// * `reader` - a mutable reference to an Iterator of Mutation. This will be consumed, it need to store Mutation in order.
/// * `snps_buffer`- a mutable reference to the VecDeque that is used as a buffer for SNPs. SNPs before the given window will
///                  be removed, the overlapping ones will be at positions 0..returned value and the first SNPs after the given window
///                  will be the last element.
pub fn find_overlapping_snps<'a, I>(window: mutations::Coordinate, reader: &mut I, snps_buffer: &mut VecDeque<&'a mutations::Mutation>) -> u32
    where I: Iterator<Item=&'a mutations::Mutation> {
    // We assume to receive ordered vcf and bed therefore we can skip vcf entries < r.start
    // if vcf > r.start+r.end empty snps_on_seq and return ---> not empty! leave them there 
    // We could use a VcfReader as reader but to be able to write unit tests more easily it is an Iterator of Mutation.
    let mut overlapping_snps = 0u32;
    let mut i = 0;
    let mut n_to_be_removed = 0;
    while i < (*snps_buffer).len() {
        let next_mut = (*snps_buffer).get(i).unwrap();
        match window.relative_position(&(*next_mut).pos) {
            mutations::Position::After => {
                n_to_be_removed += 1;
            },
            mutations::Position::Overlapping => {
                overlapping_snps += 1;
            },
            mutations::Position::Before => {
                return overlapping_snps;
            }
        }
        i += 1;        
    }
    // I am not removing them inside the previous loop to avoid borrowing issues, it is less
    // efficient though. Right now I could use a for instead of the while but will leave it 
    // in order to do the push_back inside the loop.
    while n_to_be_removed > 0 {
        let _ = (*snps_buffer).pop_front();
        n_to_be_removed -= 1;
    }
    while let Some(next_mut) = (*reader).next() {
        match window.relative_position(& next_mut.pos) {
            mutations::Position::Overlapping => {
                overlapping_snps += 1;
                (*snps_buffer).push_back(next_mut);
            },
            mutations::Position::After => {},
            mutations::Position::Before => {
                (*snps_buffer).push_back(next_mut);
                break;
            }
        }
    }
    // do we need to explicitly manage the case where VcfReader is exausted? Don't think so.
    overlapping_snps
}

pub struct MutatedSequences<'a> {
    pub genotypes: Vec<(usize, usize)>,
    pub sequences: Vec<&'a [u8]>
} 

// Probably sequences of MutatedSequences should not be a vec of slices but a vec of vecs.
pub fn obtain_seq<'a>(window: mutations::Coordinate, snps_buffer: & VecDeque<&'a mutations::Mutation>, n_overlapping: u32, reference: &'a fasta::Fasta) -> MutatedSequences<'a> {
    // snps_buffer will be empty or contain snps found in the previous window
    // if there are no overlapping snps we get the reference sequence and return only it
    // otherwise we need to build the sequences and the individual vectors.
    let ref_seq : &[u8];
    let s = window.start as usize;
    let mut e = window.end as usize;
    if e > reference.sequence.len() {
        e = reference.sequence.len();
    }
    ref_seq = &reference.sequence[s..e];
    if n_overlapping == 0 {
        MutatedSequences{ genotypes : Vec::new(), sequences: vec!(ref_seq)}
    } else {
        let n_seq = 2usize.pow(n_overlapping);
        let mut res : Vec<&'a [u8]> = Vec::with_capacity(n_seq); //Vec<&'a [u8]>
        let n_samples = 0usize; //TODO
        let genotypes : Vec<(usize, usize)> = encode_genotypes(&snps_buffer, n_overlapping, n_samples);
        let mut i = 0;
        let mut pos = s;
        let pl = vec!();
        res.push(&pl); // placeholder
        while i < n_overlapping as usize { //but what will be the order?
            let snp = snps_buffer.get(i).unwrap();
            let before_snp_seq = &reference.sequence[pos..snp.pos.start as usize];
            let mut new_res = Vec::with_capacity(n_seq); //Vec<&'a [u8]>
            for seq in res.into_iter() { 
                let mut new_seq_1 = seq.to_vec();
                let mut new_seq_2 = seq.to_vec();
                new_seq_1.extend_from_slice(before_snp_seq);
                new_seq_2.extend_from_slice(before_snp_seq);
                new_seq_1.push(snp.sequence_ref[0]);
                new_seq_2.push(snp.sequence_alt[0]);
                new_res.push(new_seq_1.as_slice());
                new_res.push(new_seq_2.as_slice());
            }
            res = new_res;
            // But is it the only thing that we can do? Ideally we have te sequences
            // already allocated inside mutations and should just point to them, but how?
            pos = snp.pos.end as usize;
            i += 1;    
        }
        // Use genotype indexes before to produce only needed sequences or?
        MutatedSequences{ genotypes : Vec::new(), sequences: res}
    }
}

pub fn encode_genotypes(snps_buffer: & VecDeque<& mutations::Mutation>, n_overlapping: u32, n_samples: usize) -> Vec<(usize, usize)> {
    //    pub genotypes: Vec<(bool, bool)>
    let mut chr_one = Vec::<String>::with_capacity(n_samples);
    let mut chr_two = Vec::<String>::with_capacity(n_samples);
    for _i in 0 .. n_samples {
        chr_one.push(String::with_capacity(n_overlapping as usize));
        chr_two.push(String::with_capacity(n_overlapping as usize));
    } 
    let mut i_snps = (n_overlapping - 1) as isize;
    while i_snps >= 0 {
        //for snp in snps_buffer.iter().rev() {
        let snp = snps_buffer.get(i_snps as usize).unwrap();
        let mut i = (n_samples - 1) as isize;
        for allele in snp.genotypes.iter() {
            match allele.0 {
                true => chr_one.get_mut(i as usize).unwrap().push('1'),
                false => chr_one.get_mut(i as usize).unwrap().push('0')
            }
            match allele.1 {
                true => chr_two.get_mut(i as usize).unwrap().push('1'),
                false => chr_two.get_mut(i as usize).unwrap().push('0')
            }
            i -= 1;
        }
        i_snps -= 1;
    }
    chr_one.iter().zip(chr_two.iter()).map(|x| (u32::from_str_radix(x.0, 2).unwrap() as usize, u32::from_str_radix(x.1, 2).unwrap() as usize)).collect()
}

// Needed structs:
// return type of get scores with bed info, n of snps, [len of individuals seqs], scores for both alleles for individuals and individual ids
// vcf entry ?