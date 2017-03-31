EXE=../target/debug/vcf_rider

OBJ=single_snp

all: $(OBJ)

clean:
	rm -rf $(OBJ)

single_snp: single_snp.vcf aaaa.pwm prova.bed ref.fasta
	$(EXE) $^ $@.bedov > $@
	diff $@.cmp $@
	diff $@.bedov.cmp $@.bedov
	rm $@.bedov
