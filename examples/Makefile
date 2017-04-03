EXE=../target/debug/vcf_rider

OBJ=single_snp small_del adj_dels

all: $(OBJ)

clean:
	rm -rf $(OBJ)

single_snp: single_snp.vcf aaaa.pwm prova.bed ref.fasta
	$(EXE) $^ $@.bedov > $@
	diff $@.cmp $@
	diff $@.bedov.cmp $@.bedov
	rm $@.bedov

small_del: smalldel.vcf aaaa.pwm prova.bed refdel.fasta
	$(EXE) $^ $@.bedov > $@
	diff $@.cmp $@
	diff $@.bedov.cmp $@.bedov
	rm $@.bedov

adj_dels: overlappingdel.vcf aaaa.pwm prova.bed ref_overlappingdel.fasta
	$(EXE) $^ $@.bedov > $@
	diff $@.cmp $@
	diff $@.bedov.cmp $@.bedov
	rm $@.bedov

