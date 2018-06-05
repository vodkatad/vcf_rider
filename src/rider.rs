use bio::io::bed;
use std::fs;
use std::io::BufWriter;
use super::fasta;
use super::mutations;
use super::indel;
use std::collections::VecDeque;
use std::io::Write;
use bit_vec::BitVec;

/// Our vcf_rider main function will receive a Vec<T: CanScoreSequence>
/// and call it for every T on subsequences of the genomes of the samples
/// doing it only for each existing subsequence once.
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
    ///           The check that sequence.len() - self.get_length() >= 0 IS NOT DONE HERE, 
    ///           we always return a score.
    /// * `sequence`- the sequence that needs to be scored, encoded as [ACGTN]-[01234]
    fn get_score(&self, pos: usize, sequence: &[u8]) -> f64;

    /// Returns the length of sequence that this object can score.
    ///
    /// # Arguments
    ///
    /// * `self` - the object with trait CanScoreSequence.
    fn get_length(&self) -> usize;

    /// Returns the name of this object.
    ///
    /// # Arguments
    ///
    /// * `self` - the object with trait CanScoreSequence.
    fn get_name(&self) -> &str;
}

/// The parameters used to setup vcf_rider are a vector of objects
/// able to score a sequence and their minimum and maximum lengths.
/// In the future it will become possible to combine scores for all the subsequences 
/// not only summing but also for example averaging them, getting the minimum or the 
/// maximum, etc.
pub struct RiderParameters<'a, T : CanScoreSequence + 'a> {
    /// The minimum length of the sequence that these T objects can score. Not used right now, we get it for 
    /// possible future optimizations.
    pub min_len: usize,
    /// The maximum length of the sequence that these T objects can score
    pub max_len: usize,
    // The vector of T objects that we want do use for scoring
    pub parameters: &'a Vec<T>
    // TODO the operation to be used to manage scores http://doc.rust-lang.org/core/ops/
}

/// The single entry point of our library, right now for ease of use in bioinformatic pipelines it simply prints the results
/// on standard output. TODO: return a suitable data structure with results.
///
/// # Arguments
///
/// * `params` - the RiderParameter that will be used to score individual sequences
/// * `vcf_path` - a path to a vcf file representing mutations on a given chr for some individuals
/// * `bed_reader` - a bed::Reader representing the genomic intervals that we want to consider. They should be on a single chr,
///                  sorted and not overlapping. They should not fall outside of the given chromosome sequence
/// * `ref_path` - the path to a fasta file with the reference sequence of the given chr
/// * `associations` - an optional string, if given this will be the file name where associations between bed ids and polymorphism will be 
///                    printed - for bed with at least an overlapping SNPs we will print the ids of the overlapping SNPs (ids are starting coord,
///                    indel_length and boolean is_indel - note that insertions are represent with a negative length and deletion with a positive one).
/// * `accept_unphased` - if unphased vcf files will be accepted. Some scores (i.e. TBA) need phased genotypes cause different phases will
///                       end up in different scores, while other (i.e. GC content) are invariant with respect to phasing and for them
///                       we will obtain sequences considering even unphased vcf as phased
///
/// # Panics 
/// If the bed defines some intervals out of the given chromosome fasta.
/// If the given fasta is empty, it's a multifasta, it's not readable or it's not a fasta.
/// It it cannot open the given associations file.
// TODO: fix check of fasta with wrong format.
pub fn get_scores<T : CanScoreSequence>(params: RiderParameters<T>, vcf_path: &str, mut bed_reader: bed::Reader<fs::File>, ref_path: &str, associations: Option<String>, accept_unphased: bool) {

    // #[cfg(debug_assertions)]   attributes on non-item statements and expressions are experimental. (see issue #15701)
    //{
    //println!("I would use {}", vcf_path);
    //println!("With parameters {} {}", params.min_len, params.max_len);
    //}

    // Load the reference fasta.
    let referenceseq: fasta::Fasta = {
        if let Ok(mut reader) = fasta::FastaReader::open_path(ref_path) {
            let referenceseq = match reader.next() {
                Some(x) => x,
                None => panic!("Empty fasta?")
            };
            let other = reader.next(); // Very inefficient to read all of it, we will read two fasta before giving the error.
            match other {
                None => (),
                Some(_) => panic!("Right now you can use this only on single chr! Multifasta format is not allowed")
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

    // Load the vcf.
    if let Ok(vcf) = mutations::VcfReader::open_path(vcf_path, accept_unphased) {
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

        // We iterate over all our bed records.
        for r in bed_reader.records() {
            let record = r.ok().expect("Error reading bed record");
            // TODO add check on the bed having ordered and not overlapping records.
            let bed_window = mutations::Coordinate{chr: "".to_owned(), start: record.start(), end: record.end()};
            // We count the number of overlapping SNPs over this bed, reading the vcf meanwhile.
            // We put in our snps_buffer all the overlapping ones and the first one outside this bed record.
            let n_overlapping = find_overlapping_snps(& bed_window, &mut vcf_reader, &mut snps_buffer);
            // We print the information about overlapping snpns in the associations file, if it has been given.
            match assoc_writer {
                Some(ref mut writer) => { print_overlapping(& snps_buffer, n_overlapping as usize, writer, &record) },
                None => {}
            }
            // We need to split the individuals/alleles in groups with uniform coordinates for this bed, i.e. with the
            // same genotypes for the overlappin gindels.
            let mut indel_manager = indel::IndelRider::new(&snps_buffer, n_overlapping, n_samples);
            // We iterate over different groups, each group is made of single chromosomes of out samples with the same
            // combination of indels genotypes for this bed.
            while let Some(chr_samples) = indel_manager.next() {
                let mut groups_snps_buffer = snps_buffer.clone();
                //println!("working on group {:?}", chr_samples);
                let mut pos = record.start();
                let mut samples : Vec<usize> = Vec::new();
                // We need to obtain the samples id for this group (XXX Do in IndelRider?)
                // chr_samples has the ids of the alleles found in this group, since we index
                // alleles for individual 1 as 0 and n_samples, for individual 2 as 1 and n_samples+1
                for allele in chr_samples.iter() {
                    samples.push(allele % n_samples); 
                }
                let n_alleles = chr_samples.len();
                //println!("{} {:?}", n_alleles, samples);
                // scores will hold the scores computed for this bed/group combination.
                // Outer vector: n_pwm, i.e. number of objects implementing CanScore that we have.
                // Inner vectors: index of the allele that has this score.
                let mut scores : Vec<Vec<f64>> = vec![vec![0f64; n_alleles]; n_pwm];
                // idx_for_seq will store correspondence between our samples
                // and the built sequences. Every tuple represent a sample (indexed by usize) 
                // and if it has the current sequence (current in the final loop done on all the possible sequences TODO).
                let mut idx_for_seq  : Vec<(usize, bool)> = Vec::<(usize, bool)>::with_capacity(n_alleles);
                while pos < record.end() {  // we do not do pos + params.max_len < r.end to avoid cumbersome management for the last portion
                    let mut wend = pos+params.max_len as u64;
                    if wend > record.end() {
                        wend = record.end();
                    }   
                    let mut window = mutations::Coordinate{chr: "".to_owned(), start: pos, end: wend};

                    // Obtain the information about the mutations that we need to manage for this group
                    // and the fixed indels: a vector of their index in groups_snps_buffer and their MutationClass
                    let mut overlapping : Vec<(usize, indel::MutationClass)> = Vec::new(); 
                    // TODO or is it better to allocate it in eccess with n overlapping capacity?
                    // This will also modify the window to access the right portion of the reference genome (longer or shorter if necessary due to indels).
                    //println!("The window was {:?}", window);
                    indel_manager.get_group_info(& mut window, & mut pos, & mut groups_snps_buffer, n_overlapping, & mut overlapping); 
                    // the group we are iterating on is known cause indel_manager is an iterator of groups.
                    //println!("And became {:?} next {}", window, pos);
                    //let n_overlapping = overlapping.iter().fold(0, |acc, &x| if x.1 == MutationClass.Manage { acc + 1} else { acc });
                    //println!("for group {:?} in window {} {} n_overlapping {}", chr_samples, window.start, window.end, n_overlapping);
                    //println!("overlapping_info {:?} ", overlapping);
                    // Obtain the encoded indexes of our genotypes, genotypes has an element for each of our samples
                    // that encodes its genotype (using only the mutation that needs to be managed here: SNPs but not indels).
                    // Every SNP status is encoded by 0 (in the BitVec) if it is reference and 1 if it is
                    // mutated. The first snp in the seq is encoded by the least significant position.
                    let genotypes : Vec<BitVec> = encode_genotypes(&groups_snps_buffer, &overlapping, &chr_samples, n_samples, &samples);
                    //println!("encoded_genotypes {:?} ", genotypes);
                
                    //println!("trying to allocate {}", 2usize.pow(n_overlapping));
                    // Obtain all the existing sequences, without repetitions, for this group in this position.
                    // TODO if we have more than XXX sequences we would break, add a check?
                    //let mut seqs : Vec<(usize, Vec<u8>)> = Vec::with_capacity(2usize.pow(n_overlapping));
                    let mut seqs : Vec<(BitVec, Vec<u8>)> = Vec::with_capacity(1usize);
                    obtain_seq(& window, & groups_snps_buffer, & overlapping, & referenceseq, & genotypes, &mut seqs, bed_window.end);
                
                    for (i, s) in seqs.into_iter() {
                        // if i in indexes genotypes -> match_indexes checks if it's there and fills a vector (idx_for_seq)
                        // with the indexes of the individuals that have the sequence that we are looping on.
                        if match_indexes(i, &mut idx_for_seq, &genotypes) {
                            // we want to compute the score for all our pwm/can_score objects on this sequence.
                            for i in 0..n_pwm {
                                let p = params.parameters.get(i).unwrap();
                                //println!("pwm name {} {} {}", p.get_name(), p.get_length(), s.len());
                                if p.get_length() <= s.len() {
                                    let score = p.get_score(0usize, &s);
                                    // now we need to sum (or smt else) the scores assigning them to the right individuals.
                                    // iterate over idx_for_seq and sum the right scores.
                                    for j in idx_for_seq.iter() {
                                        // j.0 is the wanted chr / sample index
                                        if j.1 { // if this individual [where for individuals we mean a single chr, paternal or maternal, of each sample],
                                                 // j.0, has this seq
                                                 //println!("scoring pwm {} for {}", params.parameters.get(i).unwrap().get_name(), j.0);
                                            scores[i][j.0] += score; // i indexes the pwm, j.0 the individual, two chrs are encoded by different ids.
                                        }
                                    }
                                }
                            }
                        }
                        idx_for_seq.clear();
                    }
                    //pos += 1; // if there are ins inside this window we need to avoid losing their adjacent bases.
                    // pos is managed inside get_group_info
                }

                // We print out the results for each bed when we finished processing it.
                for i in 0..n_pwm {
                    for (j, chr_sample) in chr_samples.iter().enumerate() {
                        //println!("j {} chr_s {} n_samples {}", j, chr_sample, n_samples);
                        let mut phased_allele = "allele1".to_owned();
                        if *chr_sample >= n_samples {
                            phased_allele = "allele2".to_owned();
                        }
                        let ref sample = vcf_reader.samples[(*chr_sample % n_samples)];
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

/// Function that populates a Vec<(usize, bool)> for a sequence, represented by the given index, telling us for which
/// individuals/alleles (represented by an usize) have that sequence. Returns true if at least an individual/allele
/// has this sequence.
///
/// # Arguments
///
/// * `index` - the index of the sequence that we want to consider
/// * `idx` - an empty vector to be filled with correspondences
/// * `genotypes`- a vector of the BitVec indexes of the individuals/alleles that we are studying
pub fn match_indexes(index: BitVec, idx: &mut Vec<(usize, bool)>, genotypes : &Vec<BitVec>) -> bool {
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
/// in snps_buffer all the overlapping snps and their number and then the first not overlapping snp. It uses SNPs in snps_buffer to manage
/// overlapping (always sorted on their starting coord!) bed entries. b3_e < b2_e TODO is managed/test?
///
/// # Arguments
///
/// * `window` - the window where we search for overlapping SNPs. This function should be called giving them in order with 
///              respect to genomic coords
/// * `reader` - a mutable reference to an Iterator of Mutation. This will be consumed, it needs to store Mutation objects in order
/// * `snps_buffer`- a mutable reference to the VecDeque that is used as a buffer for SNPs. SNPs before the given window will
///                  be removed, the overlapping ones will be at positions 0..returned value and the first SNPs after the given window
///                  will be the last element
pub fn find_overlapping_snps<I>(window: & mutations::Coordinate, reader: &mut I, snps_buffer: &mut VecDeque<mutations::Mutation>) -> usize
    where I: Iterator<Item=mutations::Mutation> {
    // We assume to receive ordered vcf and bed therefore we can skip vcf entries < r.start
    // if vcf > r.start+r.end empty snps_on_seq and return ---> not empty! leave them there
    // We could use a VcfReader as reader but to be able to write unit tests more easily it is an Iterator of Mutation.
    let mut overlapping_snps = 0usize;
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
                // Deletions have positive lengths and we set their end as start+1 
                // cause we do not want to get overlaps for large deletions over multiple windows
                // after this step. Our deletion will overlap a single window and be completely managed there.
                if next_mut.indel_len > 0   {
                    next_mut.pos.end = next_mut.pos.start+1;
                }
                snps_buffer.push_back(next_mut);
            },
            mutations::Position::After => {},
            mutations::Position::Before => {
                // Deletions have positive lengths and we set their end as start+1 
                // cause we do not want to get overlaps for large deletions over multiple windows
                // after this step. Our deletion will overlap a single window and be completely managed there.
                // TODO big indels starting before our bed are managed FIXME
                // TODO possible bug for deletions covering two nearby beds?
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


/// Function that prints information about the overlapping SNPs in the given writer. This is only an output function that does not compute anything
///
/// # Arguments
///
/// * `snps_buffer`- a mutable reference to the VecDeque that is used as a buffer for SNPs. We will prints all SNPs listed there save the last one
/// * `n_overlapping` - the number of SNPs in snps_buffer that overlap with our bed_record
/// * `writer` - where information about overlaps will be printed, tab delimited, as: bed name, bed start, bed end, snp_id
///    (snp_ids are starting coord, indel_length and boolean is_indel - note that insertions are represent with a negative length and deletion with a positive one)
/// * `bed_record` - the bed record whose overlaps will be printed
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

/// Function that populates a vector of sequences, given a genomic window, a buffer of SNPs overlapping with it alongside their indel status information
/// for the samples groups, the reference sequence, the genotypes encoded for our individuals and the end of the bed that we are considering.
/// The filled vector is a vector of tuples of BitVec, representing the sequences indexes, and Vec<u8>, representing the sequences itself. 
/// At this stage we will remove duplicated sequences, storing only once the ones that are the same in different individual/alleles.

/// # Arguments
/// UPTO
/// * `snps_buffer`- a mutable reference to the VecDeque that is used as a buffer for SNPs. We will prints all SNPs listed there save the last one
/// * `n_overlapping` - the number of SNPs in snps_buffer that overlap with our bed_record
/// * `writer` - where information about overlaps will be printed, tab delimited, as: bed name, bed start, bed end, snp_id
///    (snp_ids are starting coord, indel_length and boolean is_indel - note that insertions are represent with a negative length and deletion with a positive one)
/// * `bed_record` - the bed record whose overlaps will be printed
pub fn obtain_seq(window: & mutations::Coordinate, snps_buffer: & VecDeque<mutations::Mutation>, overlapping_info: & Vec<(usize, indel::MutationClass)>,
                  reference: & fasta::Fasta, genotypes : &Vec<BitVec>, seqs : &mut Vec<(BitVec, Vec<u8>)>, bed_end : u64) {
    // if there are no overlapping snps we get the reference sequence and return only it
    // otherwise we need to build the sequences
    let ref_seq : &[u8];
    let s = window.start as usize;
    let mut e = window.end as usize;
    if e as u64 > bed_end {
        e = bed_end as usize;
    }
    if e > reference.sequence.len() {
        e = reference.sequence.len();
    }
    ref_seq = &reference.sequence[s..e];
    seqs.push((BitVec::from_elem(overlapping_info.len(), false), ref_seq.to_owned()));
    //println!("non mutated window.start {} seq {:?}", s, ref_seq);
    let mut regenotypes : Vec<BitVec> = genotypes.to_vec();
    // Here we sort and dedup to avoid scoring multiple time the same sequence, we will use BitVec to make
    // samples and sequences correspond so we do not need to care about seqs vector order.
    &regenotypes.sort();
    &regenotypes.dedup();
    for encoded_geno in regenotypes.iter() {
        // We skip the all false genotype that we do not need to encode twice (already added as the reference sequence in line 368);
        if encoded_geno.none() {
            continue;
        }
        let mut seq_to_mutate = ref_seq.to_owned();
        let mut pos_adjust : isize = 0;
        for (j, &(mut_idx, ref manage)) in overlapping_info.iter().enumerate() {
            // If this polymorphism is TRUE it is mutated (==alternate allele) for this genotype
            // otherwise we have a polymorphism which is == to the reference and we do not change the sequence.
            if encoded_geno.get(j).unwrap() {
            //if ((*i) >> j) & 1 == 1 {
                let this_mut = snps_buffer.get(mut_idx as usize).unwrap();
                //println!("i {:?} j {} {:?}", encoded_geno, j, this_mut);
                match *manage {
                    // SNPs with alternate allele
                    // TODO break motivate and write: if this SNP was overlapping but...
                    indel::MutationClass::Manage(pos) => {  
                                                    let apos: usize = (pos as isize + pos_adjust) as usize; 
                                                    if apos < seq_to_mutate.len() {
                                                        seq_to_mutate[apos as usize] = this_mut.sequence_alt[0]; 
                                                    } else {
                                                        break;
                                                    }
                                                    },
                    // Insertions with alternate allele
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
                    // Deletions with alternate allele
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
                    // Polymorphism that are always == to the reference for this group
                    indel::MutationClass::Reference => { 
                                                panic!("I found smt annotated as Reference that seems mutated to me! {:?} {} i {:?} j {}", this_mut, mut_idx, encoded_geno, j);
                                                }
                }
            } 
            // The else branch here, with reference indels for this group, do not need management here for the length: this
            // is done by IndelRider. TODO CHECK.
        }
        //println!("encoded {}  window.start {} seq {:?}", i, s, seq_to_mutate);
        seqs.push((encoded_geno.clone(), seq_to_mutate));
    }
}

pub fn encode_genotypes(snps_buffer: & VecDeque<mutations::Mutation>, overlapping_info: & Vec<(usize, indel::MutationClass)>, group: &Vec<usize>, n_samples: usize, id_samples: & Vec<usize>) -> Vec<BitVec> {
    //let mut chrs : Vec<usize> = vec![0; group.len()];
    let mut chrs : Vec<BitVec> = vec![BitVec::from_elem(overlapping_info.len(), false); group.len()]; // from_elem is unstable RFC509?
    let mut bit_index = 0;
    for &(i_snp, _) in overlapping_info.iter() {
        let snp = snps_buffer.get(i_snp).unwrap();
        for (i, &i_allele) in group.iter().enumerate() {
            // chrs[i] = chrs[i] << 1;
            // we track the bit that we are modifying with bit_index and no more with shifting.
            let allele = snp.genotypes[id_samples[i]]; 
            //println!("for allele {} sample {} we see {:?}", i, id_samples[i], allele);
            if i_allele < n_samples { // snps on the first chr
                match allele.0 {
                    //true => chrs[i] |= 1,
                    true => chrs[i].set(bit_index, true),
                    false => ()
                }
            } else {  // snps on the second chr
                match allele.1 {
                    //true => chrs[i] |= 1,
                    true => chrs[i].set(bit_index, true),
                    false => ()
                }
            }
        }
        bit_index += 1;
    }
    chrs
}