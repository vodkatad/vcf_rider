EXE=../target/debug/vcf_rider

OBJ=single_snp small_del small_ins adj_dels overlapping_bed bigdel_before

all: $(OBJ)

clean:
	rm -rf $(OBJ)

single_snp: single_snp.vcf aaaa.pwm prova.bed ref.fasta bg
	$(EXE) -v $(word 1, $^) -p $(word 2, $^) -b $(word 3, $^) -r $(word 4, $^) -a $@.bedov -f $(word 5, $^) > $@
	diff $@.cmp $@
	diff $@.bedov.cmp $@.bedov
	rm $@.bedov

small_del: smalldel.vcf aaaa.pwm prova.bed refdel.fasta bg
	$(EXE) -v $(word 1, $^) -p $(word 2, $^) -b $(word 3, $^) -r $(word 4, $^) -a $@.bedov -f $(word 5, $^) > $@
	diff $@.cmp $@
	diff $@.bedov.cmp $@.bedov
	rm $@.bedov

small_ins: smallins.vcf aaaa.pwm prova.bed refins.fasta bg
	$(EXE) -v $(word 1, $^) -p $(word 2, $^) -b $(word 3, $^) -r $(word 4, $^) -a $@.bedov -f $(word 5, $^) > $@
	diff $@.cmp $@
	diff $@.bedov.cmp $@.bedov
	rm $@.bedov

adj_dels: overlappingdel.vcf aaaa.pwm prova.bed ref_overlappingdel.fasta bg
	$(EXE) -v $(word 1, $^) -p $(word 2, $^) -b $(word 3, $^) -r $(word 4, $^) -a $@.bedov -f $(word 5, $^) > $@
	diff $@.cmp $@
	diff $@.bedov.cmp $@.bedov
	rm $@.bedov

overlapping_bed: overlapping.vcf aaaa.pwm overlapping.bed ref.fasta bg
	$(EXE) -v $(word 1, $^) -p $(word 2, $^) -b $(word 3, $^) -r $(word 4, $^) -a $@.bedov -f $(word 5, $^) > $@
	diff $@.cmp $@
	diff $@.bedov.cmp $@.bedov
	rm $@.bedov

bigdel_before: bigdelbefore.vcf aaaa.pwm bigdelbefore.bed ref.fasta bg
	$(EXE) -v $(word 1, $^) -p $(word 2, $^) -b $(word 3, $^) -r $(word 4, $^) -a $@.bedov -f $(word 5, $^) > $@
	diff $@.cmp $@
	diff $@.bedov.cmp $@.bedov
	rm $@.bedov
