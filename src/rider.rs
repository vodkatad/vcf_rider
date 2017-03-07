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
        // initialize snps_buffer  VecDeque<mutations::Mutation>
        let mut snps_buffer : VecDeque<mutations::Mutation> = VecDeque::new();
        // will probably end being a VecDeque of mutations and not ref to them
        // chr check
        for r in bed_reader.records() {
            let record = r.ok().expect("Error reading record");
            // chr check
            let mut groups : Vec<u32> =  vec![0; n_samples*2]; // Vec::with_capacity(n_samples*2);
            let bed_window = mutations::Coordinate{chr: "".to_owned(), start: record.start(), end: record.end()};
            let n_overlapping = find_overlapping_snps_outer(& bed_window, &mut vcf_reader, &mut snps_buffer);      
            count_groups(&snps_buffer, n_overlapping, &mut groups, n_samples);
            let n_groups = groups.iter().max().unwrap();
            let n = *n_groups as usize;
            let mut rev_groups : Vec<Vec<u32>> = vec![Vec::new(); n+1]; // functional way to do this?
            // groups has chr samples as indexes and group ids as elements, we need to invert this array.
            for (sample, group) in groups.iter().enumerate() {
                rev_groups[*group as usize].push(sample as u32); // Mh, use all usize and stop? XXX
            }
            for chr_samples in rev_groups.iter() {
                let mut pos = record.start();
                let mut samples : Vec<u32> = Vec::new();
                for allele in chr_samples.iter() {
                    samples.push(allele % n_samples as u32);  // These are the samples id.
                }
                let n_wanted_samples = samples.len() / 2;
                //let ordered_samples = samples.iter().collect().sort();
                let n_alleles = chr_samples.len();
                //println!("{}", n_wanted_samples);
                // We change the indexing of individuals separating chrM and chrP to handle groups 
                // with different indels combo later on.
                let mut scores : Vec<Vec<f64>> = vec![vec![0f64; n_alleles]; n_pwm];
                let mut idx_for_seq  : Vec<(usize, bool)> = Vec::<(usize, bool)>::with_capacity(n_alleles);
                //let mut scores = vec![vec![(0f64); n_wanted_samples*2]; n_pwm]; // Horrible.
                while pos < record.end() {  // we do not do pos + params.max_len < r.end to avoid cumbersome management for the last portion
                    let mut wend = pos+params.max_len as u64;
                    if wend > record.end() {
                        wend = record.end();
                    }
                    let window = mutations::Coordinate{chr: "".to_owned(), start: pos, end: wend};
                    let (start_ov, end_ov) = find_overlapping_snps_inner(& window, & snps_buffer);
                    let n_overlapping = end_ov - start_ov;
                    println!("for group {:?} in window {} n_overlapping {} wend{}", chr_samples, pos,  n_overlapping, wend);
                    // TODO: pass fixed-size vector to be filled with indices.
                    let genotypes : Vec<(usize)> = encode_genotypes(&snps_buffer, start_ov, end_ov, n_wanted_samples, &samples);
                    // BEWARE: indexes of samples in the encoded genotypes are not == as the final ones.
                    let mut seqs : Vec<Vec<u8>> = Vec::with_capacity(2usize.pow(n_overlapping));
                    let len_seq = obtain_seq(& window, & snps_buffer, start_ov, end_ov, n_overlapping, & referenceseq, & genotypes, &mut seqs);
                    // TODO manage coords and len.
                    println!("for group {:?} in window {} the window is of len {}", chr_samples, pos, len_seq);
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
                        //println!("j {} chr_s {} n_wanted_samples {}", j, chr_sample, n_wanted_samples);
                        let ref sample = vcf_reader.samples[(*chr_sample % n_wanted_samples as u32) as usize];
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
pub fn find_overlapping_snps_outer<I>(window: & mutations::Coordinate, reader: &mut I, snps_buffer: &mut VecDeque<mutations::Mutation>) -> u32
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

pub fn find_overlapping_snps_inner(window: & mutations::Coordinate, snps_buffer: & VecDeque<mutations::Mutation>) -> (u32, u32) {
    let mut first_ov = 0u32;
    let mut last_ov = 0u32;
    let mut not_seen = true;
    let mut_iter : Vec<& mutations::Mutation> = snps_buffer.into_iter().collect();
    let mut i = 0;
    while let Some(next_mut) = mut_iter.get(i) {
        match window.relative_position(& next_mut.pos) {
            mutations::Position::Overlapping => {
                println!("for window {} {} seen ov {} {}", window.start, window.end, next_mut.pos.start, next_mut.pos.end);
                if not_seen {
                    first_ov = i as u32;
                    not_seen = false;
                }
                last_ov = i as u32;
            },
            mutations::Position::After => {},
            mutations::Position::Before => {
                break;
            }
        }
        i += 1;
    }
    // the second element of the tuple is the index of the first out of this window.
    if ! not_seen {
        (first_ov, last_ov+1)
    } else {
        (0, 0)
    }
}


#[allow(unused_variables)]
pub fn obtain_seq(window: & mutations::Coordinate, snps_buffer: & VecDeque<mutations::Mutation>, start_ov: u32, end_ov: u32, n_overlapping: u32,
                  reference: & fasta::Fasta, genotypes : &Vec<usize>, seqs : &mut Vec<Vec<u8>>) -> u64 {
    // snps_buffer will be empty or contain snps found in the previous window
    // if there are no overlapping snps we get the reference sequence and return only it
    // otherwise we need to build the sequences
    let ref_seq : &[u8];
    let s = window.start as usize;
    let mut e = window.end as usize;
    let mut len = window.end - window.start;
    if e > reference.sequence.len() {
        e = reference.sequence.len();
    }
    ref_seq = &reference.sequence[s..e];
    seqs.push(ref_seq.to_owned());
    //let indexes = genotypes.into_iter().somehowgetallvalues
    for i in 1..2usize.pow(n_overlapping) {
        // if i in indexes
        let mut seq_to_mutate = ref_seq.to_owned();
        let mut this_window_start = s;
        for j in start_ov .. end_ov {
            // j does not start from 0 therefore this if is not working
            if (i >> (j-start_ov)) & 1 == 1 {
                let this_mut = snps_buffer.get(j as usize).unwrap();
                if this_mut.is_indel {
                    // TODO reason about multiple indels!
                    println!("encoding {} mut.start {}  window.start {} fixed window start {} seq {:?}", i, this_mut.pos.start, s, this_window_start,  seq_to_mutate);
                    let ref mut after_mut = seq_to_mutate.split_off(this_mut.pos.start as usize - this_window_start);
                    // better to use remove/insert?
                    if this_mut.sequence_alt == vec![6u8, 6u8, 6u8] { 
                        // <DEL>, large deletion.
                        let mut del_len = this_mut.indel_len as u64;
                        //after_mut.reverse().truncate(len);
                        if del_len > after_mut.len() as u64 {
                            del_len = after_mut.len() as u64;
                        }
                        len -= del_len;
                        for k in 0 .. del_len as usize {
                            after_mut.remove(0);
                            this_window_start += 1;
                        }
                    } else { 
                        // small IN or DEL.
                        // This code should be moved inside Mutation
                        let indel_len = this_mut.indel_len as i64;
                        if indel_len > 0 {
                            // DEL
                            let mut del_len = indel_len as u64;
                            //after_mut.reverse().truncate(len);
                            if del_len > after_mut.len() as u64 {
                                del_len = after_mut.len() as u64;
                            }
                            len -= del_len;
                            for k in 0 .. del_len as usize {
                               after_mut.remove(0);
                               this_window_start += 1;
                            }
                        } else {
                            // IN
                            let in_len = -indel_len as u64;
                            len += in_len;
                            this_window_start -= in_len as usize; // overflow if bed starts at coord 0 and we have an indel there XXX TODO
                            seq_to_mutate.append(& mut this_mut.sequence_alt.to_owned());
                        }
                    }
                    seq_to_mutate.append(after_mut);
                } else {
                    seq_to_mutate[this_mut.pos.start as usize - this_window_start] = this_mut.sequence_alt[0];  
                }
                // If sequence_alt and sequence_ref will work for indels we will probably need to fill the else branch and instead of modifying the slice in place
                // adding to a String.
                // Should we worry about coords? The set of sequences returned will have different lengths but this is already managed using s.len() correctly.
                // Coords management should keep in consideration indels for determining overlapping SNPs on subsequent windows for this group in this bed.
                // Maybe this function could return the n. of nucleotide added/removed form the reference. The code to manage overlaps should consider lengths.
                // 0 works only for single SNPs, like everything else right now. // FIXME_INDELS
            }
        }
        seqs.push(seq_to_mutate);
    }
    len
}

pub fn encode_genotypes(snps_buffer: & VecDeque<mutations::Mutation>, start_ov: u32, end_ov: u32, n_samples: usize, id_samples: & Vec<u32>) -> Vec<usize> {
    let mut chrs : Vec<usize> = vec![0; n_samples*2];
    // for snp in snps_buffer.iter().rev() { // but we need to use only n_overlapping snps! 
    for i_snp in (start_ov .. end_ov).rev() {
        let snp = snps_buffer.get(i_snp as usize).unwrap();
        for i in 0 .. n_samples {
            chrs[i] = chrs[i] << 1;
            // We start indexing all M chrs starting from 0 and P from n_samples.
            // The alternative is doing i*2/i*2+1, is it more confortable? Possibly less efficient.
            chrs[i+n_samples] = chrs[i+n_samples] << 1;
            let allele = snp.genotypes[id_samples[i] as usize]; // this i is wrong, since it is indexed on the whole vcf and not our group XXX
            match allele.0 {
                true => chrs[i] |= 1,
                false => ()
            }
            match allele.1 {
                true => chrs[i+n_samples] |= 1,
                false => ()
            }
        }
    }
    chrs
}

/// Function that assigns chr samples to different groups depending on their overlapping indel alleles.
/// chr in the same group have the same alleles of the same indels, i.e. their coords are in sync.
///
/// # Arguments
///
/// * `snps_buffer`- a mutable reference to the VecDeque that is used as a buffer for SNPs. Should contain the SNPs
///    overlapping the bed that we are interested in.
/// * `n_overlapping` - the number of overlapping SNPs, since buffer will have one more.
/// * `groups` - a mutable reference to a Vec of u32 that will be filled with groups info. Indexes: samples id. Elements: groups id.
/// * `n_sample` - the number of samples (each with two alleles for each SNP) for which we have genotypes.
///
fn count_groups(snps_buffer: & VecDeque<mutations::Mutation>, n_overlapping: u32, groups: &mut Vec<u32>, n_samples: usize) {
    for (i_snp, snp) in snps_buffer.iter().enumerate() {
        if snp.is_indel && i_snp < n_overlapping as usize { // i >= n_overlapping we have finished the overlapping snps (the last one is just waiting in the buffer)
            if snp.genotypes.iter().any(|x| x.0 || x.1) {
                // we have a bisection
                //n_groups = n_groups * 2;
                // even groups have no indels at this run.
                for i_sample in 0 .. n_samples {
                    let allele = snp.genotypes[i_sample];
                    let old_group_0 = groups[i_sample];
                    let old_group_1 = groups[i_sample+n_samples];
                    groups[i_sample] = match allele.0 {
                        true => match old_group_0 {
                            0 => 1,
                            _ => 2*old_group_0+1
                        },
                        false => match old_group_0 {
                            0 => 0,
                            _ => 2*old_group_0
                        }
                    };
                    groups[i_sample+n_samples] = match allele.1 {
                        true => match old_group_1 {
                            0 => 1,
                            _ => 2*old_group_1+1,
                        },
                        false => match old_group_1 {
                            0 => 0,
                            _ => 2*old_group_1
                        }
                    };
                }
            }
        }
    }
}

// Needed structs:
// return type of get scores with bed info, n of snps, [len of individuals seqs], scores for both alleles for individuals and individual ids
// vcf entry ?

//s.set_mutation(010011001);
//s[0] -> base at index 0
