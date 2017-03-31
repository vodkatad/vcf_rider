EXE=../target/debug/vcf_rider

OBJ=single_snp

all: $(OBJ)

clean:
	rm -rf $(OBJ)
	rm *bedov

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
