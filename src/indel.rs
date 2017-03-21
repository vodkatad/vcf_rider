use super::mutations;
use std::collections::VecDeque;

// Used to classify indels and snps in groups: SNPs will be tagged as "Manage" while
// indels as "Alternative" or "Reference" according to different groups genotypes.
#[derive(Eq, PartialEq, Debug)]
pub enum MutationClass {
    Manage(usize),
    Ins(String, usize),
    Del(u32, usize)
}

pub struct IndelRider {
    groups: Vec<u32>, // groups has chr samples as indexes and group ids as elements
}

impl Iterator for IndelRider {
    type Item = Vec<u32>;

    fn next(&mut self) -> Option<Vec<u32>> {
        return None;
    }
}

impl IndelRider {
    pub fn new(snps_buffer: & VecDeque<mutations::Mutation>, n_overlapping: u32, n_samples: usize) -> IndelRider {
        let mut my_groups : Vec<u32> =  vec![0; n_samples*2];
        IndelRider::count_groups(snps_buffer, n_overlapping, & mut my_groups, n_samples);
        /*let n_groups = groups.iter().max().unwrap();
          let n = *n_groups as usize;
          let mut rev_groups : Vec<Vec<u32>> = vec![Vec::new(); n+1]; // functional way to do this?
          // groups has chr samples as indexes and group ids as elements, we need to invert this array.
          for (sample, group) in groups.iter().enumerate() {
              rev_groups[*group as usize].push(sample as u32); // Mh, use all usize and stop? XXX
          }
        */
        IndelRider{ groups: my_groups }
    }
        
    /// Function that assigns chr samples to different groups depending on their overlapping indel alleles.
    /// chr in the same group have the same alleles of the same indels, i.e. their coords are in sync.
    ///
    /// # Arguments
    ///
    /// * `snps_buffer`- a mutable reference to the VecDeque that is used as a buffer for SNPs. Should contain the SNPs
    ///    overlapping the bed that we are interested in.
    /// * `n_overlapping` - the number of overlapping SNPs, since buffer will have one more.
    /// * `groups` - a mutable reference to a Vec of u32 that will be filled with groups info. Indexes: samples id. Elements: groups id.
    /// * `n_sample` - the number of samples (each with two alleles for each SNP) for which we have genotypes.
    ///
    /// This needs also to define manage/do not manage
    fn count_groups(snps_buffer: & VecDeque<mutations::Mutation>, n_overlapping: u32, groups: &mut Vec<u32>, n_samples: usize) {
        for (i_snp, snp) in snps_buffer.iter().enumerate() {
            if snp.is_indel && i_snp < n_overlapping as usize { // i >= n_overlapping we have finished the overlapping snps (the last one is just waiting in the buffer)
                if snp.genotypes.iter().any(|x| x.0 || x.1) {
                    // we have a bisection
                    //n_groups = n_groups * 2;
                    // even groups have no indels at this run.
                    for i_sample in 0 .. n_samples {
                        let allele = snp.genotypes[i_sample];
                        let old_group_0 = groups[i_sample];
                        let old_group_1 = groups[i_sample+n_samples];
                        groups[i_sample] = match allele.0 {
                            true => match old_group_0 {
                                0 => 1,
                                _ => 2*old_group_0+1
                            },
                            false => match old_group_0 {
                                0 => 0,
                                _ => 2*old_group_0
                            }
                        };
                        groups[i_sample+n_samples] = match allele.1 {
                            true => match old_group_1 {
                                0 => 1,
                                _ => 2*old_group_1+1,
                            },
                            false => match old_group_1 {
                                0 => 0,
                                _ => 2*old_group_1
                            }
                        };
                    }
                }
            }
        }
    }
        
    // Function that given a group id and a window will return info on the overlapping SNPs for that group and on the resulting window length
    pub fn get_group_info(&self, window: & mutations::Coordinate, snps_buffer: & VecDeque<mutations::Mutation>, info: & Vec<(usize, MutationClass)>) {

    }
}