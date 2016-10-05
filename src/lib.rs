extern crate bio;


#[cfg(test)]
mod tests {
    use bio::io::bed;

    const BED_FILE: &'static [u8] = b"1\t5\t5000\tname1\tup
2\t3\t5005\tname2\tup
";

    #[test]
    fn test_reader() {
        /*let chroms = ["1", "2"];
        let starts = [5, 3];
        let ends = [5000, 5005];
        let names = ["name1", "name2"];
        let scores = ["up", "up"];*/

        let mut reader = bed::Reader::new(BED_FILE);
        for r in reader.records() {
            let record = r.ok().expect("Error reading record");
            println!("bed name: {}", record.name().expect("Error reading name"))
            /*assert_eq!(record.chrom(), chroms[i]);
            assert_eq!(record.start(), starts[i]);
            assert_eq!(record.end(), ends[i]);
            assert_eq!(record.name().expect("Error reading name"), names[i]);
            assert_eq!(record.score().expect("Error reading score"), scores[i]);*/
        }
    }
}
