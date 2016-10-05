extern crate vcf_rider;

use vcf_rider::rider;

fn main() {
    let john = rider::Person::new("John");
    john.hello();
}