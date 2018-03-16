// We need:
// Struct representing the miRNA (sequence, name and length)
// a miRNA reader
// implement CanScoreSequence for miRNAs.

// tab delimited file:
//data@tungsteno:/mnt/red/elly/bioinfotree/prj/roar/dataset/0.5/EUR/miRseeds_alteration$ head miRseeds_list 
//AAAAACC miR-548ac/548bb-3p/548d-3p/548h-3p/548z -1      hsa-miR-548ac
//AAAAACC miR-548ac/548bb-3p/548d-3p/548h-3p/548z -1      hsa-miR-548bb-3p
//AAAAACC miR-548ac/548bb-3p/548d-3p/548h-3p/548z -1      hsa-miR-548d-3p
//AAAAACC miR-548ac/548bb-3p/548d-3p/548h-3p/548z -1      hsa-miR-548h-3p
//AAAAACC miR-548ac/548bb-3p/548d-3p/548h-3p/548z -1      hsa-miR-548z
//AAAAACU miR-548aq-3p/548am-3p/548aj-3p/548ah-3p/548ae-3p/548j-3p/548x-3p        -1      hsa-miR-548ae-3p

// every subsequence score will be 1, 0:
// 1: has a seed of type 8mer, 7mer-A1 e 7mer-m8
// 0: doesn't have it
// 6mer are not considered as hits

// we always consider subsequence of length 8
// $motif = reverse $seed;
// $motif =~ tr/ACGT/TGCA/;
// $motif_6mer = substr $motif, 1, 6;
// $motif_m8 = substr $motif, 0, 1;
// $motif_A1 = 'A';

// $tentative_6mer_ref = substr $fasta_ref[$i], $j+1, 6;
// $tentative_m8_ref = substr $fasta_ref[$i], $j, 1;
// $tentative_A1_ref = substr $fasta_ref[$i], $j+7, 1;

/*        if ($tentative_6mer eq $motif_6mer)
        {
                $site = '6mer';
                $start = $j+1;
                $end = $j+6;
                if ($tentative_A1 eq $motif_A1)
                {
                        $site = '7mer-A1';
                        $end = $j+7;
                }
                if ($tentative_m8 eq $motif_m8 && $site eq '6mer')
                {
                        $site = '7mer-m8';
                        $start = $j;
                }
                if ($tentative_m8 eq $motif_m8 && $site eq '7mer-A1')
                {
                        $site = '8mer';
                        $start = $j;
                }
        }*/

/*
I get that I would need to check for == of substrings in this way:
- 7mer-A1: seq 1-6 (0based exclusive end) == to seed 6 last bases, seq 6==A
- 7mer-m8: seq 1-6 (0based exclusive end) == to seed 6 last bases, seq 0==first base of seed --> 0-7 ==seed
- 8mer: seq 1-6 (0based exclusive end) == to seed, seq 6==A, seq 0==first base of seed

So:
if seq == S1,S2,S3,S4,S5,S6,S7,A: 8mer
elseif seq start == S1,S2,S3,S4,S5,S6,S7: 7mer-m8
elseif seq end == S2,S3,S4,S5,S6,S7,A: 7mer-A1

better to compare whole string or first check the 6mer and then only single positions like the perl?

8mer if S1-S7A == my seq
7mer-A1 if S2-S7A == last 7 of my seq and first of my seq !=S1
7mer-m8 if S1-S7 == first 7 of my seed

Need to consider strand! Only the given one!
*/