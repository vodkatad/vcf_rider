# vcf_rider
 Library to efficiently compute score on individual genomes starting from vcf files 

Documentation (ongoing!) here:
https://vodkatad.github.io/vcf_rider/vcf_rider/index.html

# Installation instructions

## Install Rust (stable toolchain)

$ curl https://sh.rustup.rs -sSf | sh
(choose option 1 to install the stable toolchain by default)

Configure PATH variables in your shell without needing a new login:
$ $HOME/.cargo/env

## Clone vcf_rider repository and build all binaries

`$ clone https://github.com/vodkatad/vcf_rider.git`
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