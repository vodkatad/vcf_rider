# vcf_rider
 Library to efficiently compute score on individual genomes starting from vcf files 

Documentation (ongoing!) here:
https://vodkatad.github.io/vcf_rider/vcf_rider/index.html

## To run tests
Small number of unit tests:
`$ cargo test --lib`

More extensive tests that compare outputs with some manually checked:
`cd examples; make; make clean;`

## To build doc in target/doc
`$ cargo doc --no-deps`