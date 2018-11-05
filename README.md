# vcf_rider
 Library to efficiently compute score on individual genomes starting from vcf files 

Documentation (ongoing!) here:
https://vodkatad.github.io/vcf_rider/vcf_rider/index.html

# Installation instructions

## Install Rust (stable toolchain)

`$ curl https://sh.rustup.rs -sSf | sh`

(choose option 1 to install the stable toolchain by default)

Configure PATH variables in your shell without needing a new login:

`$ $HOME/.cargo/env`

## Clone vcf_rider repository and build all binaries

`$ git clone https://github.com/vodkatad/vcf_rider.git`

`$ cd vcf_rider`

Build the package in release mode:

`$ cargo build --release`

Binaries than can be found in ./target/release/:
`vcf_rider`
`indel_stats`
`gc_rider`
`gc_counter`

## Test your installation

Small number of unit tests:

`$ cargo test --lib`

More extensive tests that compare outputs with some manually checked:

`cargo build; cd examples; make; make clean;`

## To build doc in target/doc
`$ cargo doc --no-deps`

## Build indel stats to check if you have regions with too many mutated indels

`$ cargo build --release --bin indel_stats`

Then:

`$ ./target/release/indel_stats  your_vcf.vcf your_bed.bed`

will result in a tab delimited file with information about the overlap between regions
and mutations of the given bed and vcf, respectively. The sixth column will have true for all the regions
that right now are not correctly managed by vcf_rider, due to a huge number of indels, and should be removed
from your analyses cause they will be given incorrect scores. Future releases will fix this limitation.
