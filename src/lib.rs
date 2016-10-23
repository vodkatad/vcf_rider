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
    fn test_overlapping_inside() {
        let c1 = mutations::Coordinate{chr: "".to_owned(), start : 10, end : 15};
        let c2 = mutations::Coordinate{chr: "".to_owned(), start : 12, end : 13};
        assert_eq!(c1.relative_position(&c2), mutations::Position::Overlapping);
        assert_eq!(c2.relative_position(&c1), mutations::Position::Overlapping);
    } 
  
}
