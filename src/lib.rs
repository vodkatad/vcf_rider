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
    fn test_find_overlapping_snps() {
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
}

