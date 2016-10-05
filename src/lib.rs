extern crate bio;

pub mod rider;

#[cfg(test)]
mod tests {
    use bio::io::bed;

    // constant named BED_FILE of type pointer to vec of u8 with a 'static lifetime.
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

        // reader needs to be mut because TODO.
        let mut reader = bed::Reader::new(BED_FILE);
        for r in reader.records() {
            // .ok().expect() explanation:
            // Option<&str> is returned by record.name(). We could have used a match (TODO).
            // ok works on results, and records is an Iterator of Result<Record>
            // expect unwraps an option, yielding the content of a Some, the given &str is
            // the custom panic message if a None is there instead.

            // We rewrite the original lines as matches.
            // let record = r.ok().expect("Error reading record");
            let record = match r {
                Ok(x) => x,
                Err(e) => panic!("Error reading record {}", e),
            };
            //println!("bed name: {}", record.name().expect("Error reading name"))
            
            //let name = match record.name() {
            //    Some(x) => x,
            //    None => panic!("Error reading name"),
            //};
            //println!("bed name: {}", name)

            if let Some(name) = record.name() {
                println!("bed name: {}", name);
            } else {
                panic!("Error reading name")
            }
            /*assert_eq!(record.chrom(), chroms[i]);
            assert_eq!(record.start(), starts[i]);
            assert_eq!(record.end(), ends[i]);
            assert_eq!(record.name().expect("Error reading name"), names[i]);
            assert_eq!(record.score().expect("Error reading score"), scores[i]);*/
        }
    }
}
