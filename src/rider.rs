use bio::io::bed;
use std::fs;
use super::fasta;
use super::mutations;
use super::indel;
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
    if let Ok(vcf) = mutations::VcfReader::open_path(vcf_path, false) {
        let mut vcf_reader = vcf;
        /*for sample in & vcf_reader.samples {
            println!("sample {}", sample);
        }
        for snp in vcf {
            println!("snp {:?} {:?} {:?} {} {:?}", snp.pos, snp.sequence_ref, snp.sequence_alt, snp.id, snp.genotypes);
        }*/
        
        let n_samples = vcf_reader.samples.len();
        let n_pwm = params.parameters.len();
        let mut snps_buffer : VecDeque<mutations::Mutation> = VecDeque::new();
        for r in bed_reader.records() {
            let record = r.ok().expect("Error reading record");
            let bed_window = mutations::Coordinate{chr: "".to_owned(), start: record.start(), end: record.end()};
            let n_overlapping = find_overlapping_snps(& bed_window, &mut vcf_reader, &mut snps_buffer);
            let mut indel_manager = indel::IndelRider::new(&snps_buffer, n_overlapping, n_samples);
            // We iterate over different groups, each group is made of single chromosomes of out samples with the same
            // combination of indels genotypes for this bed.
            while let Some(chr_samples) = indel_manager.next() {
                println!("working on group {:?}", chr_samples);
                let mut pos = record.start();
                let mut samples : Vec<u32> = Vec::new();
                // We need to obtain are the samples id for this group (XXX Do in IndelRider?)
                for allele in chr_samples.iter() {
                    samples.push(allele % n_samples as u32); 
                }
                let n_alleles = chr_samples.len();
                println!("{}", n_alleles);
                // scores will hold the scores computed for this bed/group combination.
                let mut scores : Vec<Vec<f64>> = vec![vec![0f64; n_alleles]; n_pwm];
                // idx_for_seq will store correspondence between our samples
                // and the built sequences. Every tuple represent a sample (indexed by usize) 
                // and if it has the current sequence (current in the final loop done on all the possible sequences).
                // XXX TODO is this complex structure really needed?
                let mut idx_for_seq  : Vec<(usize, bool)> = Vec::<(usize, bool)>::with_capacity(n_alleles);
                while pos < record.end() {  // we do not do pos + params.max_len < r.end to avoid cumbersome management for the last portion
                    let mut wend = pos+params.max_len as u64;
                    if wend > record.end() {
                        wend = record.end();
                    }
                    let mut window = mutations::Coordinate{chr: "".to_owned(), start: pos, end: wend};

                    // Obtain the information about the mutations that we need to manage for this group
                    // and the fixed indels.
                    let mut overlapping : Vec<(usize, indel::MutationClass)> = Vec::new(); 
                    // or is it better to allocate it in eccess with n overlapping capacity?
                    // This will also modify the window to access the right portion of the reference genome (longer or shorter if necessary due to indels).
                    println!("The window was {:?}", window);
                    indel_manager.get_group_info(& mut window, &snps_buffer, n_overlapping, & mut overlapping); // He should know the group cause it is iterating on them itself.
                    println!("And became {:?}", window);
                    //let n_overlapping = overlapping.iter().fold(0, |acc, &x| if x.1 == MutationClass.Manage { acc + 1} else { acc });
                    let n_overlapping = overlapping.len() as u32;
                    println!("for group {:?} in window {} n_overlapping {} ", chr_samples, pos, n_overlapping);
                    println!("overlapping_info {:?} ", overlapping);
                    // Obtain the encoded indexes of our genotypes, genotypes has an element for each of our samples
                    // that encodes its genotype (using only the mutation that needs to be managed here, i.e. SNPs).
                    let genotypes : Vec<usize> = encode_genotypes(&snps_buffer, &overlapping, n_alleles, &samples);
                    // Obtain all the possible sequences for this group in this position.
                    let mut seqs : Vec<Vec<u8>> = Vec::with_capacity(2usize.pow(n_overlapping));
                    obtain_seq(& window, & snps_buffer, & overlapping, & referenceseq, & genotypes, &mut seqs);

                    // TODO needs updating
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
                                        // j.0 is the wanted chr / sample index
                                        if j.1 { // if this individual, j.0, has this seq
                                            scores[i][j.0] += score; // i indexes the pwm, j.0 the individual, two chrs are encoded by different ids.
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
                    //for (j, sample) in vcf_reader.samples.iter().enumerate() {
                    for (j, chr_sample) in chr_samples.iter().enumerate() {
                        println!("j {} chr_s {} n_samples {}", j, chr_sample, n_samples);
                        let ref sample = vcf_reader.samples[(*chr_sample % n_samples as u32) as usize];
                        println!("{}\t{}\t{}\t{}", record.name().expect("Error reading name"), params.parameters.get(i).unwrap().get_name(), sample, scores[i][j]);
                    }
                }
            }
        }
    }
}

pub fn match_indexes(index: usize, idx: &mut Vec<(usize, bool)>, genotypes : &Vec<usize>) -> bool {
    let mut res = false;
    for (i, allelic_tuple) in genotypes.iter().enumerate() { //could seek in a smarter way
        if *allelic_tuple == index {
            idx.push((i, true));
            res = true;
        } else {
            idx.push((i, false));
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
pub fn obtain_seq(window: & mutations::Coordinate, snps_buffer: & VecDeque<mutations::Mutation>, overlapping_info: & Vec<(usize, indel::MutationClass)>,
                  reference: & fasta::Fasta, genotypes : &Vec<usize>, seqs : &mut Vec<Vec<u8>>) {
    // snps_buffer will be empty or contain snps found in the previous window
    // if there are no overlapping snps we get the reference sequence and return only it
    // otherwise we need to build the sequences
    let n_overlapping = overlapping_info.len() as u32;
    let ref_seq : &[u8];
    let s = window.start as usize;
    let mut e = window.end as usize;
    let len = window.end - window.start;
    let n = 2usize.pow(n_overlapping);
    if e > reference.sequence.len() {
        e = reference.sequence.len();
    }
    ref_seq = &reference.sequence[s..e];
    seqs.push(ref_seq.to_owned());
    for i in 1..2usize.pow(n_overlapping) {
        if genotypes.iter().any(|&x| x == i) {
            let mut seq_to_mutate = ref_seq.to_owned();
            for (j, &(mut_idx, ref manage)) in overlapping_info.iter().enumerate() {
                if (i >> j) & 1 == 1 {
                    let this_mut = snps_buffer.get(mut_idx as usize).unwrap();
                    match *manage {
                        indel::MutationClass::Manage(pos) => seq_to_mutate[pos] = this_mut.sequence_alt[0],
                        indel::MutationClass::Ins(ref seq, pos) => {  
                                                    let ref mut after_mut = seq_to_mutate.split_off(pos);
                                                    let mut ins = seq.clone();
                                                    seq_to_mutate.append(& mut ins);
                                                    seq_to_mutate.append(after_mut);
                                                    },
                        indel::MutationClass::Del(length, pos) => {  
                                                    let ref mut after_mut = seq_to_mutate.split_off(pos);
                                                    //let ref deleted = after_mut.split_off(pos+length); TODO
                                                    seq_to_mutate.append(after_mut);
                                                    },
                        indel::MutationClass::Reference => { panic!("I found an indel annotated as Reference that seems mutated to me! {:?} {}", this_mut, mut_idx)}
                    }
                } // it is possible that we will need to manage also the else branch here, because reference indels could need management
                // to correctly manage window lenghts
            }
            println!("encoded {}  window.start {} seq {:?}", i, s, seq_to_mutate);
            seqs.push(seq_to_mutate);
        }
    }
}

pub fn encode_genotypes(snps_buffer: & VecDeque<mutations::Mutation>, overlapping_info: & Vec<(usize, indel::MutationClass)>, n_alleles: usize, id_samples: & Vec<u32>) -> Vec<usize> {
    let mut chrs : Vec<usize> = vec![0; n_alleles];
    for &(i_snp, _) in overlapping_info.iter().rev() {
        let snp = snps_buffer.get(i_snp as usize).unwrap();
        for i in 0 .. n_alleles {
            chrs[i] = chrs[i] << 1;
            let allele = snp.genotypes[id_samples[i] as usize]; 
            if i % 2 == 0 { // even if for snps on the first chr
                match allele.0 {
                    true => chrs[i] |= 1,
                    false => ()
                }
            } else {  // odd for the others
                match allele.1 {
                    true => chrs[i] |= 1,
                    false => ()
                }
            }
        }
    }
    chrs
}



// Needed structs:
// return type of get scores with bed info, n of snps, [len of individuals seqs], scores for both alleles for individuals and individual ids
// vcf entry ?

//s.set_mutation(010011001);
//s[0] -> base at index 0
