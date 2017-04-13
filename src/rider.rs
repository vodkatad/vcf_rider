use bio::io::bed;
use std::fs;
use std::io::BufWriter;
use super::fasta;
use super::mutations;
use super::indel;
use std::collections::VecDeque;
use std::io::Write;

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
pub fn get_scores<T : CanScoreSequence>(params: RiderParameters<T>, vcf_path: &str, mut bed_reader: bed::Reader<fs::File>, ref_path: &str, associations: Option<String>) {

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
    let mut assoc_writer : Option<BufWriter<fs::File>> = match associations {
         Some(x) => { 
                     let assoc_file = fs::File::create(&x).expect(&format!("Could not open {}", &x));
                     Some(BufWriter::new(assoc_file))
                    },
         None => None
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
            match assoc_writer {
                Some(ref mut writer) => { print_overlapping(& snps_buffer, n_overlapping as usize, writer, &record) },
                None => {}
            }
            let mut indel_manager = indel::IndelRider::new(&snps_buffer, n_overlapping, n_samples);
            // We iterate over different groups, each group is made of single chromosomes of out samples with the same
            // combination of indels genotypes for this bed.
            while let Some(chr_samples) = indel_manager.next() {
                let mut groups_snps_buffer = snps_buffer.clone();
                //println!("working on group {:?}", chr_samples);
                let mut pos = record.start();
                let mut samples : Vec<u32> = Vec::new();
                // We need to obtain the samples id for this group (XXX Do in IndelRider?)
                for allele in chr_samples.iter() {
                    samples.push(allele % n_samples as u32); 
                }
                let n_alleles = chr_samples.len();
                //println!("{} {:?}", n_alleles, samples);
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
                    //println!("The window was {:?}", window);
                    indel_manager.get_group_info(& mut window, & mut pos, & mut groups_snps_buffer, n_overlapping, & mut overlapping); // He should know the group cause it is iterating on them itself.
                    //println!("And became {:?} next {}", window, pos);
                    //let n_overlapping = overlapping.iter().fold(0, |acc, &x| if x.1 == MutationClass.Manage { acc + 1} else { acc });
                    //let n_overlapping = overlapping.len() as u32;
                    //println!("for group {:?} in window {} {} n_overlapping {} ", chr_samples, window.start, window.end, n_overlapping);
                    //println!("overlapping_info {:?} ", overlapping);
                    // Obtain the encoded indexes of our genotypes, genotypes has an element for each of our samples
                    // that encodes its genotype (using only the mutation that needs to be managed here, i.e. SNPs).
                    let genotypes : Vec<usize> = encode_genotypes(&groups_snps_buffer, &overlapping, &chr_samples, n_samples, &samples);
                    //println!("encoded_genotypes {:?} ", genotypes);
                    //println!("trying to allocate {}", 2usize.pow(n_overlapping));
                    // Obtain all the possible sequences for this group in this position.
                    //let mut seqs : Vec<(usize, Vec<u8>)> = Vec::with_capacity(2usize.pow(n_overlapping));
                    let mut seqs : Vec<(usize, Vec<u8>)> = Vec::new();
                    obtain_seq(& window, & groups_snps_buffer, & overlapping, & referenceseq, & genotypes, &mut seqs, bed_window.end);

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

                    for (i, s) in seqs.into_iter() {
                        // if i in indexes genotypes -> function that checks if it's there and fills a vector (idx_for_seq) with the indexes of the individuals that
                        if match_indexes(i, &mut idx_for_seq, &genotypes) {
                            // needs this score (i, 0|1)
                            for i in 0..n_pwm {
                                let p = params.parameters.get(i).unwrap();
                                //println!("pwm name {} {} {}", p.get_name(), p.get_length(), s.len());
                                if p.get_length() <= s.len() {
                                    let score = p.get_score(0usize, &s);
                                    // now we need to sum (or smt else) the scores assigning them to the right individuals.
                                    // iterate over idx_for_seq and sum the right scores.
                                    for j in idx_for_seq.iter() {
                                        // j.0 is the wanted chr / sample index
                                        if j.1 { // if this individual, j.0, has this seq
                                            //println!("scoring pwm {} for {}", params.parameters.get(i).unwrap().get_name(), j.0);
                                            scores[i][j.0] += score; // i indexes the pwm, j.0 the individual, two chrs are encoded by different ids.
                                        }
                                    }
                                }
                            }
                        }
                        idx_for_seq.clear();
                    }
                    //pos += 1; // if there are ins inside this window we need to avoid losing their adjacent bases. XXX
                    // pos is managed inside get_group_info
                }
                for i in 0..n_pwm {
                    for (j, chr_sample) in chr_samples.iter().enumerate() {
                        //println!("j {} chr_s {} n_samples {}", j, chr_sample, n_samples);
                        let mut phased_allele = "allele1".to_owned();
                        if *chr_sample as usize >= n_samples {
                            phased_allele = "allele2".to_owned();
                        }
                        let ref sample = vcf_reader.samples[(*chr_sample % n_samples as u32) as usize];
                        println!("{}\t{}\t{}\t{}\t{}\t{}\t{}", record.name().expect("Error reading name"), record.start(), record.end(),
                        params.parameters.get(i).unwrap().get_name(), sample, phased_allele, scores[i][j]);
                    }
                }
            }
        }
    }
    /*if assoc_writer.is_some() {
        assoc_writer.unwrap().flush().expect("Error while writing to the associations file");
    }*/
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
    while let Some(mut next_mut) = reader.next() {
        match window.relative_position(& next_mut.pos) {
            mutations::Position::Overlapping => {
                overlapping_snps += 1;
                if next_mut.indel_len > 0   {
                    next_mut.pos.end = next_mut.pos.start+1;
                }
                snps_buffer.push_back(next_mut);
            },
            mutations::Position::After => {},
            mutations::Position::Before => {
                if next_mut.indel_len > 0   {
                    next_mut.pos.end = next_mut.pos.start+1;
                }
                snps_buffer.push_back(next_mut);
                break;
            }
        }
    }
    // do we need to explicitly manage the case where VcfReader is exausted? Don't think so.
    overlapping_snps
}

pub fn print_overlapping(snps_buffer: & VecDeque<mutations::Mutation>, n_overlapping: usize, writer: & mut Write, bed_record: & bed::Record) {
    if n_overlapping != 0 {
        write!(writer, "{}\t{}\t{}\t", bed_record.name().expect("Error reading name"), bed_record.start(), bed_record.end()).expect("Error writing to the associations file!");
        for i in 0 .. n_overlapping-1 {
            let snp = snps_buffer.get(i).expect("error in SNPs buffer");
            write!(writer, "{}_{}_{},", snp.id, snp.indel_len, snp.is_indel).expect("Error writing to the associations file!");
        }
        let snp = snps_buffer.get(n_overlapping-1).expect("error in SNPs buffer");
        writeln!(writer, "{}_{}_{}", snp.id, snp.indel_len, snp.is_indel).expect("Error writing to the associations file!");
    }
}

pub fn obtain_seq(window: & mutations::Coordinate, snps_buffer: & VecDeque<mutations::Mutation>, overlapping_info: & Vec<(usize, indel::MutationClass)>,
                  reference: & fasta::Fasta, genotypes : &Vec<usize>, seqs : &mut Vec<(usize, Vec<u8>)>, bed_end : u64) {
    // snps_buffer will be empty or contain snps found in the previous window
    // if there are no overlapping snps we get the reference sequence and return only it
    // otherwise we need to build the sequences
    let n_overlapping = overlapping_info.len() as u32;
    let ref_seq : &[u8];
    let s = window.start as usize;
    let mut e = window.end as usize;
    if e as u64 > bed_end {
        e = bed_end as usize;
    }
    if e > reference.sequence.len() { //maybe border on bed?
        e = reference.sequence.len();
    }
    ref_seq = &reference.sequence[s..e];
    seqs.push((0, ref_seq.to_owned()));
    //println!("non mutated window.start {} seq {:?}", s, ref_seq);
    for i in 1..2usize.pow(n_overlapping) {
        if genotypes.iter().any(|&x| x == i) {
            let mut seq_to_mutate = ref_seq.to_owned();
            let mut pos_adjust : isize = 0;
            for (j, &(mut_idx, ref manage)) in overlapping_info.iter().enumerate() {
                if (i >> j) & 1 == 1 {
                    let this_mut = snps_buffer.get(mut_idx as usize).unwrap();
                    match *manage {
                        indel::MutationClass::Manage(pos) => {  
                                                        let apos: usize = (pos as isize + pos_adjust) as usize; 
                                                        if apos < seq_to_mutate.len() {
                                                            seq_to_mutate[apos as usize] = this_mut.sequence_alt[0]; 
                                                        } else {
                                                            break;
                                                        }
                                                        },
                        indel::MutationClass::Ins(ref seq, pos) => {  
                                                    //println!("managing insertion {:?}",seq_to_mutate);
                                                    let apos: usize = (pos as isize + pos_adjust) as usize; 
                                                    if apos <= seq_to_mutate.len() {
                                                        let ref mut after_mut = seq_to_mutate.split_off(apos as usize);
                                                        //println!("managing insertion {:?} {:?} {} {:?}", after_mut, seq_to_mutate, apos, seq);
                                                        let mut ins = seq.clone();
                                                        pos_adjust += ins.len() as isize;
                                                        seq_to_mutate.append(& mut ins);
                                                        seq_to_mutate.append(after_mut);
                                                    } else {
                                                        break;
                                                    }
                                                    },
                        indel::MutationClass::Del(length, pos) => {
                                                    //println!("managing del {:?}",seq_to_mutate);
                                                    let apos: usize = (pos as isize + pos_adjust) as usize; 
                                                    if apos <= seq_to_mutate.len() {
                                                        let ref mut after_mut = seq_to_mutate.split_off(apos as usize);
                                                        //println!("managing del {:?} {:?} {} {}", after_mut, seq_to_mutate, apos, length);
                                                        if length < after_mut.len() as u64 {
                                                            let ref mut after_deleted = after_mut.split_off(length as usize);
                                                            pos_adjust -= length as isize; 
                                                            seq_to_mutate.append(after_deleted);
                                                        }
                                                    } else {
                                                        break;
                                                    }
                                                    },
                       indel::MutationClass::Reference => { 
                                                    panic!("I found smt annotated as Reference that seems mutated to me! {:?} {} i {}", this_mut, mut_idx, i);
                                                    }
                    }
                } // it is possible that we will need to manage also the else branch here, because reference indels could need management
                // to correctly manage window lenghts: done by the IndelRider?
            }
            //println!("encoded {}  window.start {} seq {:?}", i, s, seq_to_mutate);
            seqs.push((i, seq_to_mutate));
        }
    }
}

pub fn encode_genotypes(snps_buffer: & VecDeque<mutations::Mutation>, overlapping_info: & Vec<(usize, indel::MutationClass)>, group: &Vec<u32>, n_samples: usize, id_samples: & Vec<u32>) -> Vec<usize> {
    let mut chrs : Vec<usize> = vec![0; group.len()];
    for &(i_snp, _) in overlapping_info.iter().rev() {
        let snp = snps_buffer.get(i_snp as usize).unwrap();
        for (i, &i_allele) in group.iter().enumerate() {
            chrs[i] = chrs[i] << 1;
            let allele = snp.genotypes[id_samples[i] as usize]; 
            //println!("for allele {} sample {} we see {:?}", i, id_samples[i], allele);
            if i_allele < n_samples as u32{ // snps on the first chr
                match allele.0 {
                    true => chrs[i] |= 1,
                    false => ()
                }
            } else {  // snps on the second chr
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
