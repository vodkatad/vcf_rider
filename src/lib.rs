//! vcf_rider: a library to efficiently compute score on individual genomes starting from vcf files 
//!
//! The idea behind vcf_rider is to exploit the fact that polymorphisms are rare. If one needs to compute
//! a sequence based score (i.e. CpG/CG content, number of PWM hits, number of miRNA seeds, Total Binding Affinity) that can
//! be computed on defined and independent windows it is not needed to reconstruct every individual genome and then compute the scores
//! separatedly. vcf_rider is able to compute scores on windows that are the same for all individuals only once. Even for polimorphic subsequences
//! it computes the scores only for the number of extant sequences and correctly assigns them to different individuals.
//! Scores for different windows can be put together in different ways - right now they are only summed but extending the lib to 
//! perform different actions should be easy via a new configuration in ther RiderParameters struct.

extern crate bio;
extern crate rust_htslib;
extern crate itertools;
extern crate bit_vec;

/// Module needed to correctly manage indels when computing scores. Indels represent a problem because they make genomes of different individuals
/// 'out of phase' and force us to divide them in different groups.
mod indel;
/// The module used for reading the fasta file representing the genome of interest.
/// Right now it should contain a single chromosome to be used with `vcf_rider`.
/// The id of the fasta should be the same used in the vcf file and with genomic regions represented in the used bed.
pub mod fasta;
/// Module representing Positional Weight Matrixes and that is able to compute their score on a given sequence.
pub mod pwm;
/// Main module of vcf_rider, its function `get_scores` is the entry point of the whole analysis.
pub mod rider;
/// Module able to load mutations from a vcf file.
pub mod mutations; 
pub mod mirna;

#[cfg(test)]
mod tests {
    use mutations::Coordinate;
    use mutations::Mutation;
    use mutations::Position;
    use indel;
    //use fasta;
    use rider;
    use std::collections::VecDeque;
    use std::iter::FromIterator;

    #[test]
    fn test_relative_position() {
        // Completely inside.
        let c1 = Coordinate{chr: "".to_owned(), start : 10, end : 15};
        let c2 = Coordinate{chr: "".to_owned(), start : 12, end : 13};
        assert_eq!(c1.relative_position(&c2), Position::Overlapping);
        assert_eq!(c2.relative_position(&c1), Position::Overlapping);
        // One sided overlap.
        let c1 = Coordinate{chr: "".to_owned(), start : 10, end : 15};
        let c2 = Coordinate{chr: "".to_owned(), start : 12, end : 20};
        assert_eq!(c1.relative_position(&c2), Position::Overlapping);
        assert_eq!(c2.relative_position(&c1), Position::Overlapping);
        // Before.
        let c1 = Coordinate{chr: "".to_owned(), start : 10, end : 15};
        let c2 = Coordinate{chr: "".to_owned(), start : 15, end : 20};
        assert_eq!(c1.relative_position(&c2), Position::Before);
        assert_eq!(c2.relative_position(&c1), Position::After);
        // ... and after (which is the same but without the end==start).
        let c1 = Coordinate{chr: "".to_owned(), start : 23, end : 24};
        let c2 = Coordinate{chr: "".to_owned(), start : 15, end : 20};
        assert_eq!(c1.relative_position(&c2), Position::After);
        assert_eq!(c2.relative_position(&c1), Position::Before);
        // Identical.
        let c1 = Coordinate{chr: "".to_owned(), start : 23, end : 24};
        let c2 = Coordinate{chr: "".to_owned(), start : 23, end : 24};
        assert_eq!(c1.relative_position(&c2), Position::Overlapping);
        assert_eq!(c2.relative_position(&c1), Position::Overlapping);
        // Other things.
        let c1 = Coordinate{chr: "".to_owned(), start : 10, end : 15};
        let c2 = Coordinate{chr: "".to_owned(), start : 9, end : 16};
        assert_eq!(c1.relative_position(&c2), Position::Overlapping);
        assert_eq!(c2.relative_position(&c1), Position::Overlapping);
    } 

       
    #[test]
    fn test_find_overlapping_snps_emptybuffer1() {
        let csnp1 = Coordinate{chr: "".to_owned(), start : 10, end : 11};
        let csnp2 = Coordinate{chr: "".to_owned(), start : 20, end : 21};
        let overlapping_mut = Mutation { id: "2".to_owned(), pos: csnp2, sequence_ref: vec!(), sequence_alt: vec!(), genotypes : vec!((true, true)), is_indel : false, indel_len: 0};
        let muts = vec!(Mutation { id: "1".to_owned(), pos: csnp1, sequence_ref: vec!(), sequence_alt: vec!(), genotypes : vec!((true, true)), is_indel : false, indel_len: 0},
                    overlapping_mut);
        let window = Coordinate{chr: "".to_owned(), start : 15, end : 30};
        let ref mut buffer = VecDeque::<Mutation>::new();
        let ref mut muts_iter = muts.into_iter();
        let n_ov = rider::find_overlapping_snps(&window, muts_iter, buffer);
        assert_eq!(n_ov, 1);
        assert_eq!(buffer.front().unwrap().id, "2".to_owned());
        assert_eq!(muts_iter.len(), 0);
    }
        
    #[test]
    fn test_find_overlapping_snps_emptybuffer2() {
        let csnp1 = Coordinate{chr: "".to_owned(), start : 10, end : 11};
        let csnp2 = Coordinate{chr: "".to_owned(), start : 20, end : 21};
        let muts = vec!(Mutation { id: "1".to_owned(), pos: csnp1, sequence_ref: vec!(), sequence_alt: vec!(), genotypes : vec!((true, true)), is_indel : false, indel_len: 0},
                    Mutation { id: "2".to_owned(), pos: csnp2, sequence_ref: vec!(), sequence_alt: vec!(), genotypes : vec!((true, true)), is_indel : false, indel_len: 0});
        let window = Coordinate{chr: "".to_owned(), start : 9, end : 30};
        let ref mut buffer = VecDeque::<Mutation>::new();
        let ref mut muts_iter = muts.into_iter();
        let n_ov = rider::find_overlapping_snps(&window, muts_iter, buffer);
        assert_eq!(n_ov, 2);
        assert_eq!(buffer.front().unwrap().id, "1".to_owned());
        assert_eq!(buffer.get(1).unwrap().id, "2".to_owned());
        assert_eq!(muts_iter.len(), 0);
    }
    
    
    #[test]
    fn test_find_overlapping_snps_emptybuffer3() {
        let csnp1 = Coordinate{chr: "".to_owned(), start : 10, end : 11};
        let csnp2 = Coordinate{chr: "".to_owned(), start : 20, end : 21};
        let muts = vec!(Mutation { id: "1".to_owned(), pos: csnp1, sequence_ref: vec!(), sequence_alt: vec!(), genotypes : vec!((true, true)), is_indel : false, indel_len: 0},
                    Mutation { id: "2".to_owned(), pos: csnp2, sequence_ref: vec!(), sequence_alt: vec!(), genotypes : vec!((true, true)), is_indel : false, indel_len: 0});
        let window = Coordinate{chr: "".to_owned(), start : 5, end : 7};
        let ref mut buffer = VecDeque::<Mutation>::new();
        let ref mut muts_iter = muts.into_iter();
        let n_ov = rider::find_overlapping_snps(&window, muts_iter, buffer);
        assert_eq!(n_ov, 0);
        assert_eq!(buffer.front().unwrap().id, "1".to_owned());
        assert_eq!(muts_iter.len(), 1);
        assert_eq!(muts_iter.next().unwrap().id, "2".to_owned());
    }

    #[test]
    fn test_find_overlapping_snps_buffer_overlap() {
        let csnp1 = Coordinate{chr: "".to_owned(), start : 10, end : 11};
        let csnp2 = Coordinate{chr: "".to_owned(), start : 20, end : 21};
        let muts = Vec::new();
        let mut1 = Mutation { id: "1".to_owned(), pos: csnp1, sequence_ref: vec!(), sequence_alt: vec!(), genotypes : vec!((true, true)), is_indel : false, indel_len: 0};
        let mut2 = Mutation { id: "2".to_owned(), pos: csnp2, sequence_ref: vec!(), sequence_alt: vec!(), genotypes : vec!((true, true)), is_indel : false, indel_len: 0};
        let ref mut buffer = VecDeque::from(vec!(mut1, mut2));
        let window = Coordinate{chr: "".to_owned(), start : 5, end : 30};
        let ref mut muts_iter = muts.into_iter();
        let n_ov = rider::find_overlapping_snps(&window, muts_iter, buffer);
        assert_eq!(n_ov, 2);
        assert_eq!(buffer.front().unwrap().id, "1".to_owned());
        assert_eq!(buffer.get(1).unwrap().id, "2".to_owned());
        assert_eq!(muts_iter.len(), 0);
    }
    
    #[test]
    fn test_find_overlapping_snps_buffer_after() {
        let csnp1 = Coordinate{chr: "".to_owned(), start : 10, end : 11};
        let csnp2 = Coordinate{chr: "".to_owned(), start : 20, end : 21};
        let muts = Vec::new();
        let mut1 = Mutation { id: "1".to_owned(), pos: csnp1, sequence_ref: vec!(), sequence_alt: vec!(), genotypes : vec!((true, true)), is_indel : false, indel_len: 0};
        let mut2 = Mutation { id: "2".to_owned(), pos: csnp2, sequence_ref: vec!(), sequence_alt: vec!(), genotypes : vec!((true, true)), is_indel : false, indel_len: 0};
        let ref mut buffer = VecDeque::from(vec!(mut1, mut2));
        let window = Coordinate{chr: "".to_owned(), start : 5, end : 7};
        let ref mut muts_iter = muts.into_iter();
        let n_ov = rider::find_overlapping_snps(&window, muts_iter, buffer);
        assert_eq!(n_ov, 0);
        assert_eq!(buffer.front().unwrap().id, "1".to_owned());
        assert_eq!(buffer.get(1).unwrap().id, "2".to_owned());
        assert_eq!(muts_iter.len(), 0);
    }
    
    #[test]
    fn test_find_overlapping_snps_buffer_after_removing() {
        let csnp0 = Coordinate{chr: "".to_owned(), start : 3, end : 4};
        let csnp1 = Coordinate{chr: "".to_owned(), start : 10, end : 11};
        let csnp2 = Coordinate{chr: "".to_owned(), start : 20, end : 21};
        let muts = Vec::new();
        let mut0 = Mutation { id: "0".to_owned(), pos: csnp0, sequence_ref: vec!(), sequence_alt: vec!(), genotypes : vec!((true, true)), is_indel : false, indel_len: 0};
        let mut1 = Mutation { id: "1".to_owned(), pos: csnp1, sequence_ref: vec!(), sequence_alt: vec!(), genotypes : vec!((true, true)), is_indel : false, indel_len: 0};
        let mut2 = Mutation { id: "2".to_owned(), pos: csnp2, sequence_ref: vec!(), sequence_alt: vec!(), genotypes : vec!((true, true)), is_indel : false, indel_len: 0};
        let ref mut buffer = VecDeque::from(vec!(mut0, mut1, mut2));
        let window = Coordinate{chr: "".to_owned(), start : 5, end : 7};
        let ref mut muts_iter = muts.into_iter();
        let n_ov = rider::find_overlapping_snps(&window, muts_iter, buffer);
        assert_eq!(n_ov, 0);
        assert_eq!(buffer.front().unwrap().id, "1".to_owned());
        assert_eq!(buffer.get(1).unwrap().id, "2".to_owned());
        assert_eq!(muts_iter.len(), 0);
    }
    
    #[test]
    fn test_find_overlapping_snps_buffer_before() {
        let csnp1 = Coordinate{chr: "".to_owned(), start : 10, end : 11};
        let csnp2 = Coordinate{chr: "".to_owned(), start : 20, end : 21};
        let muts = Vec::new();
        let mut1 = Mutation { id: "1".to_owned(), pos: csnp1, sequence_ref: vec!(), sequence_alt: vec!(), genotypes : vec!((true, true)), is_indel : false, indel_len: 0};
        let mut2 = Mutation { id: "2".to_owned(), pos: csnp2, sequence_ref: vec!(), sequence_alt: vec!(), genotypes : vec!((true, true)), is_indel : false, indel_len: 0};
        let ref mut buffer = VecDeque::from(vec!(mut1, mut2));
        let window = Coordinate{chr: "".to_owned(), start : 25, end : 37};
        let ref mut muts_iter = muts.into_iter();
        let n_ov = rider::find_overlapping_snps(&window, muts_iter, buffer);
        assert_eq!(n_ov, 0);
        assert_eq!(buffer.len(), 0);
        assert_eq!(muts_iter.len(), 0);
    }
    
    #[test]
    fn test_find_overlapping_snps_buffer_mixed() {
        let csnp1 = Coordinate{chr: "".to_owned(), start : 10, end : 11};
        let csnp2 = Coordinate{chr: "".to_owned(), start : 20, end : 21};
        let csnp3 = Coordinate{chr: "".to_owned(), start : 23, end : 24};
        let mut1 = Mutation { id: "1".to_owned(), pos: csnp1, sequence_ref: vec!(), sequence_alt: vec!(), genotypes : vec!((true, true)), is_indel : false, indel_len: 0};
        let mut2 = Mutation { id: "2".to_owned(), pos: csnp2, sequence_ref: vec!(), sequence_alt: vec!(), genotypes : vec!((true, true)), is_indel : false, indel_len: 0};
        let mut buffer = VecDeque::from(vec!(mut1, mut2)); 
        let muts = vec!(Mutation { id: "3".to_owned(), pos: csnp3, sequence_ref: vec!(), sequence_alt: vec!(), genotypes : vec!((true, true)), is_indel : false, indel_len: 0});
        let mut muts_iter = muts.into_iter();
        let window = Coordinate{chr: "".to_owned(), start : 15, end : 37};
        let n_ov = rider::find_overlapping_snps(&window, &mut muts_iter, &mut buffer);
        assert_eq!(n_ov, 2);
        assert_eq!(muts_iter.len(), 0);
        assert_eq!(buffer.len(), 2);
        assert_eq!(buffer.front().unwrap().id, "2".to_owned());
        assert_eq!(buffer.get(1).unwrap().id, "3".to_owned());
    }
    // TODO: some tests calling find_overlapping_snps several times. 
    // We do not need to test it on indels, since at this level everything is the same as for SNPs.

    /*
    #[test]
    fn test_obtain_seq() {
        let window = Coordinate{chr: "".to_owned(), start : 0, end : 2};
        let ref reference = fasta::Fasta{id: "1".to_owned(), sequence: vec!(0,1,2,3), background : vec!(0.298947240099661, 0.200854143743417, 0.200941012710477, 0.299257603446445)};
        let ref mut buffer = VecDeque::<Mutation>::new();
        let mut seqs : Vec<Vec<u8>> = Vec::with_capacity(1);
        rider::obtain_seq(&window, &buffer, 0, &reference, &vec!(), &mut seqs);
        assert_eq!(seqs[0], [0, 1]);
        assert_eq!(seqs.len(), 1);
    }

    #[test]
    fn test_encode_genotypes() {
        let csnp1 = Coordinate{chr: "".to_owned(), start : 10, end : 11};
        let csnp2 = Coordinate{chr: "".to_owned(), start : 20, end : 21};
        let mut1 = Mutation { id: "1".to_owned(), pos: csnp1, sequence_ref: vec!(), sequence_alt: vec!(), genotypes : vec!((false, true)), is_indel : false};
        let mut2 = Mutation { id: "2".to_owned(), pos: csnp2, sequence_ref: vec!(), sequence_alt: vec!(), genotypes : vec!((true, false)), is_indel : false};
        let buffer = VecDeque::from(vec!(mut1, mut2));
        let indexes = rider::encode_genotypes(&buffer, 2u32, 1usize);
        // We have one individual, two overlapping snps and our guy is indexed as 10 ->  2 and 01 -> 1
        assert_eq!(indexes, vec!((2, 1)));
    }

    #[test]
    fn test_encode_genotypes_single_snp() {
        let csnp1 = Coordinate{chr: "".to_owned(), start : 10, end : 11};
        let csnp2 = Coordinate{chr: "".to_owned(), start : 20, end : 21};
        let mut1 = Mutation { id: "1".to_owned(), pos: csnp1, sequence_ref: vec!(), sequence_alt: vec!(), genotypes : vec!((false, true)), is_indel : false};
        let mut2 = Mutation { id: "2".to_owned(), pos: csnp2, sequence_ref: vec!(), sequence_alt: vec!(), genotypes : vec!((true, false)), is_indel : false};
        let buffer = VecDeque::from(vec!(mut1, mut2));
        let indexes = rider::encode_genotypes(&buffer, 1u32, 1usize);
        // We have one individual, one overlapping snps and our guy is indexed as 0 -> 0 and 1 -> 1
        assert_eq!(indexes, vec!((0, 1)));
    }

    #[test]
    fn test_obtain_seq_two_snps() {
        let csnp1 = Coordinate{chr: "".to_owned(), start : 1, end : 2};
        let csnp2 = Coordinate{chr: "".to_owned(), start : 3, end : 4};
        let mut1 = Mutation { id: "1".to_owned(), pos: csnp1, sequence_ref: vec!(1), sequence_alt: vec!(3), genotypes : vec!((false, true)), is_indel : false}; // T-C
        let mut2 = Mutation { id: "2".to_owned(), pos: csnp2, sequence_ref: vec!(2), sequence_alt: vec!(0), genotypes : vec!((true, false)), is_indel : false}; // G-A
        let buffer = VecDeque::from(vec!(mut1, mut2));
        let reference = fasta::Fasta{id: "".to_owned(), sequence: vec!(0,1,0,2), background : vec!()}; // ATAG
        let window = Coordinate{chr: "".to_owned(), start: 0, end: 4};
        let n_overlapping = 2u32;
        let mut seqs : Vec<Vec<u8>> = Vec::with_capacity(2usize.pow(n_overlapping));
        rider::obtain_seq(&window, &buffer, n_overlapping, &reference, &vec!((2,1)), &mut seqs);
        assert_eq!(seqs[0], vec!(0,1,0,2));
        assert_eq!(seqs[1], vec!(0,3,0,2));
        assert_eq!(seqs[2], vec!(0,1,0,0));
        assert_eq!(seqs[3], vec!(0,3,0,0));
        assert_eq!(seqs.len(), 4);
    }

    #[test]
    fn test_obtain_seq_two_snps_2() {
        // Just to test window positioning.
        let csnp1 = Coordinate{chr: "".to_owned(), start : 11, end : 12};
        let csnp2 = Coordinate{chr: "".to_owned(), start : 13, end : 14};
        let mut1 = Mutation { id: "1".to_owned(), pos: csnp1, sequence_ref: vec!(1), sequence_alt: vec!(3), genotypes : vec!((false, true)), is_indel : false}; // T-C
        let mut2 = Mutation { id: "2".to_owned(), pos: csnp2, sequence_ref: vec!(2), sequence_alt: vec!(0), genotypes : vec!((true, false)), is_indel : false}; // G-A
        let buffer = VecDeque::from(vec!(mut1, mut2));
        let reference = fasta::Fasta{id: "".to_owned(), sequence: vec!(0,0,0,0,0,0,0,0,0,0,0,1,0,2), background : vec!()}; // ATAG
        let window = Coordinate{chr: "".to_owned(), start: 10, end: 14};
        let n_overlapping = 2u32;
        let mut seqs : Vec<Vec<u8>> = Vec::with_capacity(2usize.pow(n_overlapping));
        rider::obtain_seq(&window, &buffer, n_overlapping, &reference, &vec!((2,1)), &mut seqs);
        assert_eq!(seqs[0], vec!(0,1,0,2));
        assert_eq!(seqs[1], vec!(0,3,0,2));
        assert_eq!(seqs[2], vec!(0,1,0,0));
        assert_eq!(seqs[3], vec!(0,3,0,0));
        assert_eq!(seqs.len(), 4);
    }*/

    #[test]
    fn test_indel_stats() {
        let csnp1 = Coordinate{chr: "".to_owned(), start : 10, end : 11};
        let csnp2 = Coordinate{chr: "".to_owned(), start : 20, end : 21};
        let csnp3 = Coordinate{chr: "".to_owned(), start : 30, end : 31};
        let muts = vec!(Mutation { id: "2".to_owned(), pos: csnp1, sequence_ref: vec!(), sequence_alt: vec!(), 
                        genotypes : vec!((true, true), (true, false), (false, false),(false, true)), is_indel : true, indel_len: 0},
                        Mutation { id: "1".to_owned(), pos: csnp2, sequence_ref: vec!(), sequence_alt: vec!(), 
                        genotypes : vec!((true, false), (false, false), (false, false),(false, false)), is_indel : true, indel_len: 0},
                        Mutation { id: "2".to_owned(), pos: csnp3, sequence_ref: vec!(), sequence_alt: vec!(), genotypes : vec!((true, true)), is_indel : false, indel_len: 0});
        let ref mut buffer = VecDeque::<Mutation>::from_iter(muts.into_iter());
        let n_samples = 4;
        let mut groups : Vec<usize> =  vec![0; n_samples*2];
        indel::IndelRider::count_groups(buffer, 2, &mut groups, n_samples);
        assert_eq!(groups, vec![3, 2, 0, 0, 2, 0, 0, 2]);
    }

    #[test]
    fn test_indel_stats_2() {
        let csnp1 = Coordinate{chr: "".to_owned(), start : 10, end : 11};
        let csnp2 = Coordinate{chr: "".to_owned(), start : 10, end : 11};
        let csnp3 = Coordinate{chr: "".to_owned(), start : 10, end : 11};
        let csnp4 = Coordinate{chr: "".to_owned(), start : 10, end : 11};
        let csnp5 = Coordinate{chr: "".to_owned(), start : 10, end : 11};
        let csnp6 = Coordinate{chr: "".to_owned(), start : 10, end : 11};
        // TOdo RECHECK with the groups tree and individuals the wanted results and genotype tuples!
        let muts = vec!(Mutation { id: "1".to_owned(), pos: csnp1, sequence_ref: vec!(), sequence_alt: vec!(), 
                        genotypes : vec!((true, false), (true, false), (true, false), (true, false)), is_indel : true, indel_len: 0},
                        Mutation { id: "2".to_owned(), pos: csnp2, sequence_ref: vec!(), sequence_alt: vec!(), 
                        genotypes : vec!((true, true), (true, true), (false, false),(false, false)), is_indel : true, indel_len: 0},
                        Mutation { id: "3".to_owned(), pos: csnp3, sequence_ref: vec!(), sequence_alt: vec!(), 
                        genotypes : vec!((false, false), (true, true), (false, false),(true, true)), is_indel : true, indel_len: 0},
                        Mutation { id: "4".to_owned(), pos: csnp4, sequence_ref: vec!(), sequence_alt: vec!(), genotypes : vec!(), is_indel : false, indel_len: 0},
                        Mutation { id: "5".to_owned(), pos: csnp5, sequence_ref: vec!(), sequence_alt: vec!(), genotypes : vec!(), is_indel : false, indel_len: 0},
                        Mutation { id: "6".to_owned(), pos: csnp6, sequence_ref: vec!(), sequence_alt: vec!(), 
                        genotypes : vec!(), is_indel : true, indel_len: 0},
                    );
        // Coords are always the same and genotypes for snps are wrong but ideally we should never look at them.
        // The lst indel should be skipped (it is the first Mutation outside our window in this setup.
        let ref mut buffer = VecDeque::<Mutation>::from_iter(muts.into_iter());
        let n_samples = 4;
        let mut groups : Vec<usize> =  vec![0; n_samples*2];
        indel::IndelRider::count_groups(buffer, 5, &mut groups, n_samples);
        assert_eq!(groups, vec![6, 7, 4, 5, 2, 3, 0, 1]); // NAY
    }
} 


