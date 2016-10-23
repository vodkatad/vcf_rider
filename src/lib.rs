extern crate bio;
extern crate rust_htslib;
extern crate itertools;

pub mod fasta;
pub mod pwm;
pub mod rider;
pub mod mutations; 

#[cfg(test)]
mod tests {
    use mutations;

    #[test]
    fn test_relative_position() {
        // Completely inside.
        let c1 = mutations::Coordinate{chr: "".to_owned(), start : 10, end : 15};
        let c2 = mutations::Coordinate{chr: "".to_owned(), start : 12, end : 13};
        assert_eq!(c1.relative_position(&c2), mutations::Position::Overlapping);
        assert_eq!(c2.relative_position(&c1), mutations::Position::Overlapping);
        // One sided overlap.
        let c1 = mutations::Coordinate{chr: "".to_owned(), start : 10, end : 15};
        let c2 = mutations::Coordinate{chr: "".to_owned(), start : 12, end : 20};
        assert_eq!(c1.relative_position(&c2), mutations::Position::Overlapping);
        assert_eq!(c2.relative_position(&c1), mutations::Position::Overlapping);
        // Before.
        let c1 = mutations::Coordinate{chr: "".to_owned(), start : 10, end : 15};
        let c2 = mutations::Coordinate{chr: "".to_owned(), start : 15, end : 20};
        assert_eq!(c1.relative_position(&c2), mutations::Position::Before);
        assert_eq!(c2.relative_position(&c1), mutations::Position::After);
        // ... and after (which is the same but without the end==start).
        let c1 = mutations::Coordinate{chr: "".to_owned(), start : 23, end : 24};
        let c2 = mutations::Coordinate{chr: "".to_owned(), start : 15, end : 20};
        assert_eq!(c1.relative_position(&c2), mutations::Position::After);
        assert_eq!(c2.relative_position(&c1), mutations::Position::Before);
        // Identical.
        let c1 = mutations::Coordinate{chr: "".to_owned(), start : 23, end : 24};
        let c2 = mutations::Coordinate{chr: "".to_owned(), start : 23, end : 24};
        assert_eq!(c1.relative_position(&c2), mutations::Position::Overlapping);
        assert_eq!(c2.relative_position(&c1), mutations::Position::Overlapping);
        // Other things.
        let c1 = mutations::Coordinate{chr: "".to_owned(), start : 10, end : 15};
        let c2 = mutations::Coordinate{chr: "".to_owned(), start : 9, end : 16};
        assert_eq!(c1.relative_position(&c2), mutations::Position::Overlapping);
        assert_eq!(c2.relative_position(&c1), mutations::Position::Overlapping);
    } 
  
}
