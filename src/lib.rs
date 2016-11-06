extern crate bio;
extern crate rust_htslib;
extern crate itertools;

pub mod fasta;
pub mod pwm;
pub mod rider;
pub mod mutations; 

#[cfg(test)]
mod tests {
    use mutations::Coordinate;
    use mutations::Mutation;
    use mutations::Position;
    use fasta;
    use rider;
    use std::collections::VecDeque;

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
        let muts = vec!(Mutation { id: "1".to_owned(), pos: csnp1, sequence_ref: vec!(), sequence_alt: vec!(), genotypes : vec!((true, true))},
                    Mutation { id: "2".to_owned(), pos: csnp2, sequence_ref: vec!(), sequence_alt: vec!(), genotypes : vec!((true, true))});
        let window = Coordinate{chr: "".to_owned(), start : 15, end : 30};
        let ref mut buffer = VecDeque::<&Mutation>::new();
        let ref mut muts_iter = muts.iter();
        let n_ov = rider::find_overlapping_snps(window, muts_iter, buffer);
        assert_eq!(n_ov, 1);
        assert_eq!(buffer.front().unwrap().id, muts[1].id);
        assert_eq!(muts_iter.len(), 0);
    }
    
    #[test]
    fn test_find_overlapping_snps_emptybuffer2() {
        let csnp1 = Coordinate{chr: "".to_owned(), start : 10, end : 11};
        let csnp2 = Coordinate{chr: "".to_owned(), start : 20, end : 21};
        let muts = vec!(Mutation { id: "1".to_owned(), pos: csnp1, sequence_ref: vec!(), sequence_alt: vec!(), genotypes : vec!((true, true))},
                    Mutation { id: "2".to_owned(), pos: csnp2, sequence_ref: vec!(), sequence_alt: vec!(), genotypes : vec!((true, true))});
        let window = Coordinate{chr: "".to_owned(), start : 9, end : 30};
        let ref mut buffer = VecDeque::<&Mutation>::new();
        let ref mut muts_iter = muts.iter();
        let n_ov = rider::find_overlapping_snps(window, muts_iter, buffer);
        assert_eq!(n_ov, 2);
        assert_eq!(buffer.front().unwrap().id, muts[0].id);
        assert_eq!(buffer.get(1).unwrap().id, muts[1].id);
        assert_eq!(muts_iter.len(), 0);
    }

    #[test]
    fn test_find_overlapping_snps_emptybuffer3() {
        let csnp1 = Coordinate{chr: "".to_owned(), start : 10, end : 11};
        let csnp2 = Coordinate{chr: "".to_owned(), start : 20, end : 21};
        let muts = vec!(Mutation { id: "1".to_owned(), pos: csnp1, sequence_ref: vec!(), sequence_alt: vec!(), genotypes : vec!((true, true))},
                    Mutation { id: "2".to_owned(), pos: csnp2, sequence_ref: vec!(), sequence_alt: vec!(), genotypes : vec!((true, true))});
        let window = Coordinate{chr: "".to_owned(), start : 5, end : 7};
        let ref mut buffer = VecDeque::<&Mutation>::new();
        let ref mut muts_iter = muts.iter();
        let n_ov = rider::find_overlapping_snps(window, muts_iter, buffer);
        assert_eq!(n_ov, 0);
        assert_eq!(buffer.front().unwrap().id, muts[0].id);
        assert_eq!(muts_iter.len(), 1);
        assert_eq!(muts_iter.next().unwrap().id, muts[1].id);
    }

    #[test]
    fn test_find_overlapping_snps_buffer_overlap() {
        let csnp1 = Coordinate{chr: "".to_owned(), start : 10, end : 11};
        let csnp2 = Coordinate{chr: "".to_owned(), start : 20, end : 21};
        let muts = Vec::new();
        let mut1 = Mutation { id: "1".to_owned(), pos: csnp1, sequence_ref: vec!(), sequence_alt: vec!(), genotypes : vec!((true, true))};
        let mut2 = Mutation { id: "2".to_owned(), pos: csnp2, sequence_ref: vec!(), sequence_alt: vec!(), genotypes : vec!((true, true))};
        let ref mut buffer = VecDeque::from(vec!(&mut1, &mut2));
        let window = Coordinate{chr: "".to_owned(), start : 5, end : 30};
        let ref mut muts_iter = muts.iter();
        let n_ov = rider::find_overlapping_snps(window, muts_iter, buffer);
        assert_eq!(n_ov, 2);
        assert_eq!(buffer.front().unwrap().id, mut1.id);
        assert_eq!(buffer.get(1).unwrap().id, mut2.id);
        assert_eq!(muts_iter.len(), 0);
    }

    #[test]
    fn test_find_overlapping_snps_buffer_after() {
        let csnp1 = Coordinate{chr: "".to_owned(), start : 10, end : 11};
        let csnp2 = Coordinate{chr: "".to_owned(), start : 20, end : 21};
        let muts = Vec::new();
        let mut1 = Mutation { id: "1".to_owned(), pos: csnp1, sequence_ref: vec!(), sequence_alt: vec!(), genotypes : vec!((true, true))};
        let mut2 = Mutation { id: "2".to_owned(), pos: csnp2, sequence_ref: vec!(), sequence_alt: vec!(), genotypes : vec!((true, true))};
        let ref mut buffer = VecDeque::from(vec!(&mut1, &mut2));
        let window = Coordinate{chr: "".to_owned(), start : 5, end : 7};
        let ref mut muts_iter = muts.iter();
        let n_ov = rider::find_overlapping_snps(window, muts_iter, buffer);
        assert_eq!(n_ov, 0);
        assert_eq!(buffer.front().unwrap().id, mut1.id);
        assert_eq!(buffer.get(1).unwrap().id, mut2.id);
        assert_eq!(muts_iter.len(), 0);
    }


    #[test]
    fn test_find_overlapping_snps_buffer_before() {
        let csnp1 = Coordinate{chr: "".to_owned(), start : 10, end : 11};
        let csnp2 = Coordinate{chr: "".to_owned(), start : 20, end : 21};
        let muts = Vec::new();
        let mut1 = Mutation { id: "1".to_owned(), pos: csnp1, sequence_ref: vec!(), sequence_alt: vec!(), genotypes : vec!((true, true))};
        let mut2 = Mutation { id: "2".to_owned(), pos: csnp2, sequence_ref: vec!(), sequence_alt: vec!(), genotypes : vec!((true, true))};
        let ref mut buffer = VecDeque::from(vec!(&mut1, &mut2));
        let window = Coordinate{chr: "".to_owned(), start : 25, end : 37};
        let ref mut muts_iter = muts.iter();
        let n_ov = rider::find_overlapping_snps(window, muts_iter, buffer);
        assert_eq!(n_ov, 0);
        assert_eq!(buffer.len(), 0);
        assert_eq!(muts_iter.len(), 0);
    }

    #[test]
    fn test_find_overlapping_snps_buffer_mixed() {
        let csnp1 = Coordinate{chr: "".to_owned(), start : 10, end : 11};
        let csnp2 = Coordinate{chr: "".to_owned(), start : 20, end : 21};
        let csnp3 = Coordinate{chr: "".to_owned(), start : 23, end : 24};
        let mut1 = Mutation { id: "1".to_owned(), pos: csnp1, sequence_ref: vec!(), sequence_alt: vec!(), genotypes : vec!((true, true))};
        let mut2 = Mutation { id: "2".to_owned(), pos: csnp2, sequence_ref: vec!(), sequence_alt: vec!(), genotypes : vec!((true, true))};
        
        //let mut buffer = VecDeque::from(vec!(&mut1, &mut2)); 
        // If it is here it does not compile
        // this is due to the lifetimes declared in find_overlapping_snps.
        let muts = vec!(Mutation { id: "3".to_owned(), pos: csnp3, sequence_ref: vec!(), sequence_alt: vec!(), genotypes : vec!((true, true))});
        let mut muts_iter = muts.iter();
        let mut buffer = VecDeque::from(vec!(&mut1, &mut2)); // Here it works.

        let window = Coordinate{chr: "".to_owned(), start : 15, end : 37};
        let n_ov = rider::find_overlapping_snps(window, &mut muts_iter, &mut buffer);
        assert_eq!(n_ov, 2);
        assert_eq!(muts_iter.len(), 0);
        assert_eq!(buffer.len(), 2);
        assert_eq!(buffer.front().unwrap().id, mut2.id);
        assert_eq!(buffer.get(1).unwrap().id, muts[0].id);
    }
    // TODO: some tests calling find_overlapping_snps several times.

    #[test]
    fn test_obtain_seq() {
        let window = Coordinate{chr: "".to_owned(), start : 0, end : 2};
        let ref reference = fasta::Fasta{id: "1".to_owned(), sequence: vec!(0,1,2,3), background : vec!(0.298947240099661, 0.200854143743417, 0.200941012710477, 0.299257603446445)};
        let ref mut buffer = VecDeque::<Mutation>::new();
        let mut seqs : Vec<Vec<u8>> = Vec::with_capacity(1);
        rider::obtain_seq(window, &buffer, 0, &reference, vec!(), &mut seqs);
        assert_eq!(seqs[0], [0, 1]);
        assert_eq!(seqs.len(), 1);
    }

    #[test]
    fn test_encode_genotypes() {
        let csnp1 = Coordinate{chr: "".to_owned(), start : 10, end : 11};
        let csnp2 = Coordinate{chr: "".to_owned(), start : 20, end : 21};
        let mut1 = Mutation { id: "1".to_owned(), pos: csnp1, sequence_ref: vec!(), sequence_alt: vec!(), genotypes : vec!((false, true))};
        let mut2 = Mutation { id: "2".to_owned(), pos: csnp2, sequence_ref: vec!(), sequence_alt: vec!(), genotypes : vec!((true, false))};
        let buffer = VecDeque::from(vec!(&mut1, &mut2));
        let indexes = rider::encode_genotypes(&buffer, 2u32, 1usize);
        // We have one individual, two overlapping snps and our guy is indexed as 10 ->  2 and 01 -> 1
        assert_eq!(indexes, vec!((2, 1)));
    }

    #[test]
    fn test_obtain_seq_two_snps() {
        let csnp1 = Coordinate{chr: "".to_owned(), start : 1, end : 2};
        let csnp2 = Coordinate{chr: "".to_owned(), start : 3, end : 4};
        let mut1 = Mutation { id: "1".to_owned(), pos: csnp1, sequence_ref: vec!(1), sequence_alt: vec!(3), genotypes : vec!((false, true))}; // T-C
        let mut2 = Mutation { id: "2".to_owned(), pos: csnp2, sequence_ref: vec!(2), sequence_alt: vec!(0), genotypes : vec!((true, false))}; // G-A
        let buffer = VecDeque::from(vec!(mut1, mut2));
        let reference = fasta::Fasta{id: "".to_owned(), sequence: vec!(0,1,0,2), background : vec!()}; // ATAG
        let window = Coordinate{chr: "".to_owned(), start: 0, end: 4};
        let n_overlapping = 2u32;
        let mut seqs : Vec<Vec<u8>> = Vec::with_capacity(2usize.pow(n_overlapping));
        rider::obtain_seq(window, &buffer, n_overlapping, &reference, vec!((2,1)), &mut seqs);
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
        let mut1 = Mutation { id: "1".to_owned(), pos: csnp1, sequence_ref: vec!(1), sequence_alt: vec!(3), genotypes : vec!((false, true))}; // T-C
        let mut2 = Mutation { id: "2".to_owned(), pos: csnp2, sequence_ref: vec!(2), sequence_alt: vec!(0), genotypes : vec!((true, false))}; // G-A
        let buffer = VecDeque::from(vec!(mut1, mut2));
        let reference = fasta::Fasta{id: "".to_owned(), sequence: vec!(0,0,0,0,0,0,0,0,0,0,0,1,0,2), background : vec!()}; // ATAG
        let window = Coordinate{chr: "".to_owned(), start: 10, end: 14};
        let n_overlapping = 2u32;
        let mut seqs : Vec<Vec<u8>> = Vec::with_capacity(2usize.pow(n_overlapping));
        rider::obtain_seq(window, &buffer, n_overlapping, &reference, vec!((2,1)), &mut seqs);
        assert_eq!(seqs[0], vec!(0,1,0,2));
        assert_eq!(seqs[1], vec!(0,3,0,2));
        assert_eq!(seqs[2], vec!(0,1,0,0));
        assert_eq!(seqs[3], vec!(0,3,0,0));
        assert_eq!(seqs.len(), 4);
    }
}


