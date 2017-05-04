    use super::mutations;
    use std::collections::VecDeque;
    use std::collections::HashMap;
    
    // Used to classify indels and snps in groups: SNPs will be tagged as "Manage" while
    // indels with Ins or Del (if this group has their alternative allele).
    // String (Vec<u8>) is the sequence that needs to be inserted for ins, usize is the coord inside the window
    // and u64 the length of deletions.
    #[derive(Eq, PartialEq, Debug)]
    pub enum MutationClass {
        Manage(usize),
        Ins(Vec<u8>, usize),
        Del(u64, usize),
        Reference
    }

    pub struct IndelRider {
        groups: HashMap<usize, Vec<usize>>, // Vec of usize for samples indexes, usize keys are groups ids.
        indexes: Vec<usize>,
        next_group: usize,
        n_samples_tot: usize
    }

    impl Iterator for IndelRider {
        type Item = Vec<usize>;

        fn next(&mut self) -> Option<Vec<usize>> {
            if self.next_group == self.indexes.len() {
                None
            } else {
                let next_group_i = self.indexes[self.next_group];
                let res = self.groups.get(& next_group_i).unwrap();
                self.next_group += 1;
                Some(res.to_vec()) //XXX FIXME but why do we need to clone it? Do we need lifetimes or...?
            }
        }
    }

    impl IndelRider {
        pub fn new(snps_buffer: & VecDeque<mutations::Mutation>, n_overlapping: u32, n_samples: usize) -> IndelRider {
            let mut groups : Vec<usize> =  vec![0; n_samples*2];
            IndelRider::count_groups(snps_buffer, n_overlapping, & mut groups, n_samples);
            let mut rev_groups : HashMap<usize, Vec<usize>> = HashMap::new();
            //println!("cgroups {:?}", groups);
            // groups has chr samples as indexes and group ids as elements, we need to invert our perspective to iterate easily on groups.
            for (sample, group) in groups.iter().enumerate() {
                //rev_groups[*group as usize].push(sample as u32); // Mh, use all usize and stop? XXX
                if rev_groups.contains_key(group) {
                    let mut samples = rev_groups.get_mut(group).unwrap();
                    samples.push(sample);
                } else {
                    rev_groups.insert(*group, vec![sample]);
                }
            }
            //println!("cgroups {:?}", rev_groups);
            let groups_indexes : Vec<usize> = rev_groups.keys().map(|&x| x).collect(); 
            IndelRider{ groups: rev_groups, indexes: groups_indexes, next_group: 0, n_samples_tot: n_samples}
            //  IndelRider{ groups: rev_groups, next_group: 0, n_samples_tot: n_samples}
        }
            
        /// Function that assigns chr samples to different groups depending on their overlapping indel alleles.
        /// chr in the same group have the same alleles of the same indels, i.e. their coords are in sync.
        /// CAUTION: if there are more than 64 indels it will incurr in overflow errors! FIXME
        ///
        /// # Arguments
        ///
        /// * `snps_buffer`- a mutable reference to the VecDeque that is used as a buffer for SNPs. Should contain the SNPs
        ///    overlapping the bed that we are interested in.
        /// * `n_overlapping` - the number of overlapping SNPs, since buffer will have one more.
        /// * `groups` - a mutable reference to a Vec of usize that will be filled with groups info. Indexes: samples id. Elements: groups id.
        /// * `n_sample` - the number of samples (each with two alleles for each SNP) for which we have genotypes.
        ///
        pub fn count_groups(snps_buffer: & VecDeque<mutations::Mutation>, n_overlapping: u32, groups: &mut Vec<usize>, n_samples: usize) {
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
        pub fn get_group_info(&self, window: & mut mutations::Coordinate, next_pos: & mut u64, snps_buffer: & mut VecDeque<mutations::Mutation>, n_overlapping: u32, info: & mut Vec<(usize, MutationClass)>) {
            // Right now the logic is a bit twisted cause we change coords for snps when we get a deletion but we change window.end for overlapping indels...
            // I got why I was changing in ends...to catch their overlap across window borders, but that is wrong. We need to find a way to manage indels across window borders.
            let mut pos_managed: bool = false;
            let mut last_del_start_win: bool = false;
            for (i_snp, snp) in snps_buffer.iter_mut().enumerate() {
                if i_snp < n_overlapping as usize { // i >= n_overlapping we have finished the overlapping snps (the last one is just waiting in the buffer)
                    let next_group_i = self.indexes[self.next_group-1];
                    let this_g_v = self.groups.get(&(next_group_i)).unwrap();
                    let mut group_genotypes : Vec<bool> = Vec::with_capacity(this_g_v.len());
                    for i_sample in this_g_v.to_vec() {
                        let index = i_sample as usize % self.n_samples_tot;
                        if i_sample < self.n_samples_tot {
                            group_genotypes.push(snp.genotypes.get(index).unwrap().0);
                            //println!("for sample {} seen {}", i_sample, snp.genotypes.get(index).unwrap().0);
                        } else {
                            group_genotypes.push(snp.genotypes.get(index).unwrap().1);
                            //println!("for sample {} seen {}", i_sample, snp.genotypes.get(index).unwrap().1);
                        }
                    }
                    let mut res_mutclass = MutationClass::Manage(0); // the majority are SNPs so we start with this.
                    let snp_coords = mutations::Coordinate{ chr: snp.pos.chr.to_owned(), start: snp.pos.start, end: snp.pos.end };
                    let mut snp_end_overlap_borders = snp.pos.end;
                    let mut len_modifier : i64 = 0;
                    //println!("group {:?} genotypes {:?}", self.groups[self.next_group-1], group_genotypes);
                    // We fix coords for snps that comes after a deletion.
                    if group_genotypes.iter().any(|&x| x) {
                        if snp.is_indel {
                            if snp.indel_len > 0 {
                                // we need to know how to move coords of SNPs after deletion of this bed: no! we modify our window end
                                // therefore we do not need to move SNPs around.
                                snp_end_overlap_borders += snp.indel_len as u64;
                            } else {
                                snp_end_overlap_borders += (-snp.indel_len) as u64;
                            }
                            //for ins we get less reference since we have inserted bases for this group (snp.len is negative for ins)
                            //for del we need to get more reference since we have removed bases.   
                            len_modifier = snp.indel_len; // we do not modify the window here 
                        }
                    } else {
                        res_mutclass = MutationClass::Reference;
                        // we have a SNP always reference in this group or an indel always reference.
                    }
                    // Determine overlap
                    // We need to use a window with a modified end that considers all indels, Before and Overlapping -> but only to define its 
                    // start, the length will be changed only considering Overlapping indels.
                    let sub_window = mutations::Coordinate{ chr: window.chr.to_owned(), start: window.start, end: window.end};
                    if snp_end_overlap_borders != snp.pos.end { // since we consider all SNPs of length 1
                        snp_end_overlap_borders -= 1;
                    }
                    let snp_coords_overlap = mutations::Coordinate{ chr: snp_coords.chr.to_owned(), start: snp_coords.start, end: snp_end_overlap_borders};
                    match snp_coords_overlap.relative_position_overlap(&sub_window) {
                        (mutations::Position::Before, _) => {   //println!("seen {} before", snp_coords.start)
                                                            },
                        (mutations::Position::Overlapping, overlap) => { 
                                                            if snp.indel_len != 0 { // else it is an exhausted insertion and it does not overlap anymore.
                                                                let ov = overlap.unwrap();
                                                                //println!("snp {} {} {} {}", snp_coords.start, snp_coords.end, window.start, window.end);
                                                                //println!("ov {} {}", ov.start, ov.end);
                                                                let pos = (ov.start-window.start) as usize;
                                                                let mut ov_len_modifier = (ov.end - ov.start) as u64;
                                                                if len_modifier < 0 {
                                                                    window.end -= ov_len_modifier as u64;
                                                                    let mut ins = snp.sequence_alt.to_owned();
                                                                    if snp_coords.start == ov.start && pos == 0 {
                                                                        // we do not have to modify the window start otherwise we risk getting wrong sequences.
                                                                        // we change our next_window and "eat out" the insertion step by step?
                                                                        // no it does not work cause it will always overlap at least a base?
                                                                        snp.indel_len += 1;
                                                                        snp.sequence_alt.remove(0);
                                                                        if snp_coords.end > ov.end && snp.indel_len == (snp_coords.end - ov.end) as i64 { 
                                                                            //manage long insertion
                                                                            snp.pos.start = ov.end;  
                                                                        } 
                                                                        pos_managed = true; 
                                                                        *next_pos = snp.pos.start;
                                                                    }
                                                                    if ov.start > snp_coords.start {
                                                                        for _r in snp_coords.start .. ov.start {
                                                                            ins.remove(0);
                                                                        }
                                                                    } else if ov.end < snp_coords_overlap.end {
                                                                        ins.split_off((ov.end - snp_coords_overlap.start) as usize);
                                                                    }
                                                                    res_mutclass = MutationClass::Ins(ins, pos);
                                                                    last_del_start_win = false;
                                                                } else if len_modifier > 0 {
                                                                    if last_del_start_win && *next_pos == snp_coords_overlap.end { 
                                                                        // very funny case with two adjacent deletions.
                                                                        pos_managed = true;
                                                                        *next_pos = snp_coords_overlap.end + 1; // are we skipping smt?
                                                                    }
                                                                    if pos == 0 { 
                                                                        // this del starts with this window, the next window should start
                                                                        // right after it.
                                                                        last_del_start_win = true;
                                                                        pos_managed = true;
                                                                        *next_pos = snp_coords_overlap.end + 1; // are we skipping smt?
                                                                    } else {
                                                                        last_del_start_win = false;
                                                                    }
                                                                    if snp_coords_overlap.end <= window.end { 
                                                                        // if it is > this del eats all this window so we do not want to enlarge it otherwise we risk getting 
                                                                        // unexistent reference bases for this deletion.
                                                                        window.end += ov_len_modifier as u64;                                                                        
                                                                    } else {
                                                                        // we enlarge it past the whole deletion but eat all of it
                                                                        window.end = snp_coords_overlap.end + ov_len_modifier; // + ov_len_modifier to keep window of the same lengths
                                                                        ov_len_modifier = snp.indel_len as u64;
                                                                    }
                                                                    // del length should be changed if they are across the window
                                                                    res_mutclass = MutationClass::Del(ov_len_modifier, pos);
                                                                } else if res_mutclass != MutationClass::Reference {
                                                                    res_mutclass = MutationClass::Manage(pos);
                                                                    last_del_start_win = false;
                                                                }
                                                                info.push((i_snp, res_mutclass));
                                                            }
                                                        },
                        (mutations::Position::After, _) => {    //println!("seen {} after", snp_coords.start); 
                                                                break } 
                    }
                }
            }
            if ! pos_managed {
                *next_pos += 1;
            }
        }
    }