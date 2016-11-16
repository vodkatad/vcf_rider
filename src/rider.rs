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
    //println!("I would use {}", vcf_path);
    //println!("With parameters {} {}", params.min_len, params.max_len);
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

    //println!("Fasta ref {}", referenceseq.id);
    // load vcf -> open file, skip # headers, first real entry
    // We could use a VcfReader similar to others.
    if let Ok(vcf) = mutations::VcfReader::open_path(vcf_path) {
        let mut vcf_reader = vcf;
        /*for sample in & vcf_reader.samples {
            println!("sample {}", sample);
        }
        for snp in vcf {
            println!("snp {:?} {:?} {:?} {} {:?}", snp.pos, snp.sequence_ref, snp.sequence_alt, snp.id, snp.genotypes);
        }*/
        
        let n_samples = vcf_reader.samples.len();
        let n_pwm = params.parameters.len();
        let mut scores : Vec<Vec<(f64, f64)>> = vec![vec![(0f64, 0f64); n_samples]; n_pwm];
        let mut idx_for_seq  : Vec<(usize, (bool,bool))> = Vec::<(usize, (bool,bool))>::with_capacity(n_samples);
        // initialize snps_buffer  VecDeque<mutations::Mutation>
        let mut snps_buffer : VecDeque<mutations::Mutation> = VecDeque::new();
        // will probably end being a VecDeque of mutations and not ref to them
        // chr check
        for r in bed_reader.records() {
            let record = r.ok().expect("Error reading record");
            //println!("bed name: {}", record.name().expect("Error reading name"));
            // chr check
            let mut pos = record.start();
            while pos < record.end() {  // we do not do pos + params.max_len < r.end to avoid cumbersome management for the last portion
		let mut wend = pos+params.max_len as u64;
		if wend > record.end() {
			wend = record.end();
		}
                let window = mutations::Coordinate{chr: "".to_owned(), start: pos, end: wend};
                let n_overlapping = find_overlapping_snps(& window, &mut vcf_reader, &mut snps_buffer);
                // TODO: pass fixed-size vector to be filled with indices.
                let genotypes : Vec<(usize, usize)> = encode_genotypes(&snps_buffer, n_overlapping, n_samples);
                let mut seqs : Vec<Vec<u8>> = Vec::with_capacity(2usize.pow(n_overlapping));
                obtain_seq(& window, & snps_buffer, n_overlapping, & referenceseq, & genotypes, &mut seqs);
                //println!("genotypes {:?}", genotypes);
                // this will give us 2^n seqs where n in the n of snps found in r.start-r.rstart+params.max_len
                // seqs will be ordered in a specific order: the first one is the reference one and the last one
                // the one with all mutated alleles.  Every SNP status is encoded by 0 if it is reference and 1 if it is
                // mutated. The first snp in the seq is encoded by the least significant position.
                // In this way it will be easy to build for each individual chr the indexes linking
                // them to their sequences. There will be a vector of tuples (usize,usize) for this.
                // obtain seq will return sequences of length params.max_len if possible otherwise shorter
                // ones and we will check to call get_score only on the right parameters
                // for every p params.parameters call on seq
                for (i, s) in seqs.iter().enumerate() {
                    //println!("seq {} {:?}", i, s);
                    // if i in indexes genotypes -> function that checks if it's there and fills a vector (idx_for_seq) with the indexes of the individuals that
                    if match_indexes(i, &mut idx_for_seq, &genotypes) {
                        // needs this score (i, 0|1)
                        for i in 0..n_pwm {
                            let p = params.parameters.get(i).unwrap();
                            //println!("pwm name {} {} {}", p.get_name(), p.get_length(), s.len());
                            if p.get_length() <= s.len() {
                                let score = p.get_score(0usize, s);
                                // now we need to sum (or smt else) the scores assigning them to the right individuals.
                                // iterate over idx_for_seq and sum the right scores.
                                for j in idx_for_seq.iter() {
                                    // j.0 is the wanted samples index, j.1.0 and 0.1 the info about the two chromosomes
                                    if (j.1).0 {
                                        scores[i][j.0].0 += score;
                                    } 
                                    if (j.1).1 {
                                        scores[i][j.0].1 += score;
                                    }
                                }
                            }
                        }
                    }
                    idx_for_seq.clear();
                }
                pos += 1;
            }
            for i in 0..n_pwm {
                for (j, sample) in vcf_reader.samples.iter().enumerate() {
                    println!("{}\t{}\t{}\t{}\t{}", record.name().expect("Error reading name"), params.parameters.get(i).unwrap().get_name(), sample, scores[i][j].0, scores[i][j].1);
                }
            }
            scores = vec![vec![(0f64, 0f64); n_samples]; n_pwm]; // Horrible.
        }
    }
}

pub fn match_indexes(index: usize, idx: &mut Vec<(usize, (bool,bool))>, genotypes : &Vec<(usize, usize)>) -> bool {
    let mut res = false;
    for (i, allelic_tuple) in genotypes.iter().enumerate() { //could seek in a smarter way
        match (allelic_tuple.0, allelic_tuple.1) {
            (i1, i2) if i1 == index && i2 != index => { idx.push((i, (true, false))); res = true;},
            (i1, i2) if i1 != index && i2 == index => { idx.push((i, (false, true))); res = true;},
            (i1, i2) if i1 == index && i2 == index => { idx.push((i, (true, true))); res = true;},
            _ => () // rething about this 
        }
    }
    return res;
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
pub fn find_overlapping_snps<I>(window: & mutations::Coordinate, reader: &mut I, snps_buffer: &mut VecDeque<mutations::Mutation>) -> u32
    where I: Iterator<Item=mutations::Mutation> {
    // We assume to receive ordered vcf and bed therefore we can skip vcf entries < r.start
    // if vcf > r.start+r.end empty snps_on_seq and return ---> not empty! leave them there
    // We could use a VcfReader as reader but to be able to write unit tests more easily it is an Iterator of Mutation.
    let mut overlapping_snps = 0u32;
    let mut i = 0;
    let mut n_to_be_removed = 0;
    let mut window_before_next_snp = false;
    while i < snps_buffer.len() {
        let next_mut = snps_buffer.get(i).unwrap();
        match window.relative_position(& next_mut.pos) {
            mutations::Position::After => {
                n_to_be_removed += 1;
            },
            mutations::Position::Overlapping => {
                overlapping_snps += 1;
            },
            mutations::Position::Before => {
                window_before_next_snp = true;
            }
        }
        i += 1;
    }
    // I am not removing them inside the previous loop to avoid borrowing issues, it is less
    // efficient though. Right now I could use a for instead of the while but will leave it
    // in order to do the push_back inside the loop.
    /*while n_to_be_removed > 0 {
        let _ = snps_buffer.pop_front(); 
        n_to_be_removed -= 1;
    }*/
    snps_buffer.drain(0..n_to_be_removed);
    // split_off does not seem ok, maybe drain?
    // test efficiency of split_off, but it does new allocations
    //let new_snp_b = snps_buffer.split_off(n_to_be_removed as usize);
    //*snps_buffer = new_snp_b;
    if window_before_next_snp {
        return overlapping_snps
    }
    while let Some(next_mut) = reader.next() {
        match window.relative_position(& next_mut.pos) {
            mutations::Position::Overlapping => {
                overlapping_snps += 1;
                snps_buffer.push_back(next_mut);
            },
            mutations::Position::After => {},
            mutations::Position::Before => {
                snps_buffer.push_back(next_mut);
                break;
            }
        }
    }
    // do we need to explicitly manage the case where VcfReader is exausted? Don't think so.
    overlapping_snps
}

#[allow(unused_variables)]
pub fn obtain_seq(window: & mutations::Coordinate, snps_buffer: & VecDeque<mutations::Mutation>, n_overlapping: u32,
                  reference: & fasta::Fasta, genotypes : &Vec<(usize, usize)>, seqs : &mut Vec<Vec<u8>>) {
    // snps_buffer will be empty or contain snps found in the previous window
    // if there are no overlapping snps we get the reference sequence and return only it
    // otherwise we need to build the sequences
    let ref_seq : &[u8];
    let s = window.start as usize;
    let mut e = window.end as usize;
    if e > reference.sequence.len() {
        e = reference.sequence.len();
    }
    ref_seq = &reference.sequence[s..e];
    seqs.push(ref_seq.to_owned());
    //let indexes = genotypes.into_iter().somehowgetallvalues
    for i in 1..2usize.pow(n_overlapping) {
        // if i in indexes
        let mut seq_to_mutate = ref_seq.to_owned();
        for j in 0 .. n_overlapping {
            if (i >> j) & 1 == 1 {
                let this_mut = snps_buffer.get(j as usize).unwrap();
                seq_to_mutate[this_mut.pos.start as usize - s] = this_mut.sequence_alt[0]; 
                // 0 works only for single SNPs, like everything else right now.
            }
        }
        seqs.push(seq_to_mutate);
    }
}

pub fn encode_genotypes(snps_buffer: & VecDeque<mutations::Mutation>, n_overlapping: u32, n_samples: usize) -> Vec<(usize, usize)> {
    let mut chr_one : Vec<usize> = vec![0; n_samples];
    let mut chr_two : Vec<usize> = vec![0; n_samples];
    // for snp in snps_buffer.iter().rev() { // but we need to use only n_overlapping snps! 
    for i_snp in (0 .. n_overlapping).rev() {
        let snp = snps_buffer.get(i_snp as usize).unwrap();
        for i in 0 .. n_samples {
            chr_one[i] = chr_one[i] << 1;
            chr_two[i] = chr_two[i] << 1;
            let allele = snp.genotypes[i];
            match allele.0 {
                true => chr_one[i] |= 1,
                false => ()
            }
            match allele.1 {
                true => chr_two[i] |= 1,
                false => ()
            }
        }
    }
    chr_one.into_iter().zip(chr_two.into_iter()).map(|x| (x.0, x.1)).collect()
}

// Needed structs:
// return type of get scores with bed info, n of snps, [len of individuals seqs], scores for both alleles for individuals and individual ids
// vcf entry ?

//s.set_mutation(010011001);
//s[0] -> base at index 0
