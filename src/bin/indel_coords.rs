extern crate vcf_rider;
extern crate bio;
extern crate itertools; 

use vcf_rider::mutations::*;
use vcf_rider::fasta;
use vcf_rider::mutations;
use std::collections::VecDeque;

fn main() {
    let csnp1 = Coordinate{chr: "".to_owned(), start : 10, end : 11};
    let csnp2 = Coordinate{chr: "".to_owned(), start : 14, end : 15};
    let csnp3 = Coordinate{chr: "".to_owned(), start : 19, end : 20};
    let csnp4 = Coordinate{chr: "".to_owned(), start : 190, end : 200};
    let mut1 = Mutation { id: "1".to_owned(), pos: csnp1, sequence_ref: vec!(0), 
                        sequence_alt: vec!(0,3,3), genotypes : vec!((false, false)), 
                        is_indel : true, indel_len : -2}; // A  ATT
    let mut2 = Mutation { id: "2".to_owned(), pos: csnp2, sequence_ref: vec!(0,0,0), 
                        sequence_alt: vec!(0), genotypes : vec!((false, false)), 
                        is_indel : true, indel_len : 2}; // AAA A
    let mut3 = Mutation { id: "3".to_owned(), pos: csnp3, sequence_ref: vec!(0), 
                        sequence_alt: vec!(6,6,6), genotypes : vec!((false, false)), 
                        is_indel : true, indel_len : 6}; // A	<DEL>, END=26
    let mut4 = Mutation { id: "4".to_owned(), pos: csnp4, sequence_ref: vec!(0), 
                        sequence_alt: vec!(6,6,6), genotypes : vec!((false, false)), 
                        is_indel : true, indel_len : 6}; // we do not care.

    let buffer = VecDeque::from(vec!(mut1, mut2, mut3, mut4));
    let reference = fasta::Fasta{id: "".to_owned(), sequence: vec![0; 200], background : vec!()};
    let window = Coordinate{chr: "".to_owned(), start: 10, end: 16};
    //let mut seqs : Vec<Vec<u8>> = Vec::with_capacity(2usize.pow(n_overlapping)); -> this needs to be created with the right length in the function.
    let mut seqs : Vec<Vec<u8>> = Vec::new();
    let end_coords = obtain_seqs(&window, &buffer, &reference, &vec!(0,0,0,0), &mut seqs);
    println!("{:?}", seqs);
    println!("{:?}", end_coords);
    // We want to obtain 2^2 seqs on this window:
    // no mutations, 00, end coord 16
    // vec![0; 6]
    // mut1 with alt allele, 01, end coord 14
    // vec![0,3,3,0,0,0]
    // 10, end coord 17 (15 and 16 are gone due to the indel)
    // vec![0,0,0,0,0,0] // we have a deleted 0, maybe change the seq to visually check that it is happened?
    // 11, end coord is 13, do not overlap with the other mut...so we have different codes...it is impossible, we will need to work KNOWING that we 
    // have a single group!
    // vec![0,3,3,0,0,0]
}

/// Black box that determines which mutations in snps_buffer need to be considered for the given window
/// and builds all the possible sequences. The resulting sequences will have the same lengths even if we find
/// indels and the function will return the coordinates on the reference that each of them end at (do we need them?).
/// Different lengths should only occur at the end of the reference seq. 
/// This also needs to put inside
/// genotypes the encoded genotypes for the overlapping mutations therefore groups management is also due?
/// Or another box? But it would need to find _again_ the overlapping SNPs.
/// 
/// NOTE: right now this will be called for a single group, therefore all the "used" sequences will have
/// the same coordinates on the reference for sequences ends. Reason about two possibilities:
/// - remove group identification and management, because it seems that having this function that does
///   the work could leave not much to do outside for coords management...or not, because the windows
///   that needs to be created will diverge in term of coords!
/// - pass info to this function about the genotypes of the group being analyzed (their indexes should do)
///   and build only the needed sequences, that will be different only for SNPs and not indels.
///   I AM CONVINCED THAT we will need to work KNOWING that we 
///   have a single group!
/// # Arguments
///
/// * `window` - the window where we search for overlapping SNPs. This function should be called giving them in order.
/// * `snps_buffer`- a mutable reference to the VecDeque that is used as a buffer for SNPs. SNPs before the given window will
///                  be removed, the overlapping ones will be at positions 0..returned value and the first SNPs after the given window
///                  will be the last element.
/// * `reference` - the reference sequence, a reference to a Fasta.
/// * `genotypes` - right now unused.
/// * `seqs` - a Vector of all the possible sequences generated by the mutations falling on the given window.
#[allow(unused_variables)]
pub fn obtain_seqs(window: & mutations::Coordinate, snps_buffer: & VecDeque<mutations::Mutation>,
                  reference: & fasta::Fasta, genotypes : &Vec<usize>, seqs : &mut Vec<Vec<u8>>) -> Vec<u64> {
    vec![0u64]

}