################################################################################
# Construct raw FSTs from Rscripts. TODO: implement in c++                     #
# Construct fst/coati.fst                                                      #
# Original code: Reed A. Cartwright https://github.com/reedacartwright/toycoati#
################################################################################

RSCRIPT = Rscript --vanilla
SYMS = fst/dna_syms.txt
FASTA = fasta/

# Construct an FST for a MG94 codon model
fst/mutation.fst: scripts/mutation.R $(SYMS)
	@$(RSCRIPT) $< | fstcompile --arc_type=standard --isymbols=$(SYMS) \
		--osymbols=$(SYMS) - | fstrmepsilon | fstarcsort --sort_type=olabel > $@

# Construct an FST for a geometric indel model
fst/indel.fst: scripts/indel.R $(SYMS)
	@$(RSCRIPT) $< | fstcompile --arc_type=standard --isymbols=$(SYMS) \
		--osymbols=$(SYMS) - | fstrmepsilon | fstarcsort --sort_type=ilabel > $@

################################################################################
# Marginalized mutation model
################################################################################

# Construct an FST for a codon-marginal MG94 codon model
fst/marg_pos.fst: scripts/marg_mutation.R $(SYMS)
	@$(RSCRIPT) $< 1 | fstcompile --arc_type=standard --isymbols=$(SYMS) \
		--osymbols=$(SYMS) - | fstrmepsilon | fstarcsort --sort_type=ilabel > $@

# Construct an FST for a nucleotide-marginal MG94 codon model
fst/dna_marg.fst: scripts/marg_mutation.R $(SYMS)
	@$(RSCRIPT) $< 2 | fstcompile --arc_type=standard --isymbols=$(SYMS) \
		--osymbols=$(SYMS) - | fstrmepsilon | fstarcsort --sort_type=ilabel > $@

################################################################################
# Find alignment graph & its shortest path
################################################################################

aln/%: build/coati  fst/indel.fst fst/marg_pos.fst fst/dna_marg.fst fst/mutation.fst
	@echo "Aligning "$*
	@./build/coati -f fasta/$* -m toy-marginal -o aln/$* #-w aln/weights.csv

################################################################################
# Drawing and printing FST                                                     #
################################################################################

# create a graph of an FST
%.dot: %.fst
	@fstdraw --isymbols=$(SYMS) --osymbols=$(SYMS) $< > $@

# draw an FST
%.pdf: %.dot
	@dot -Tpdf -o $@ $<

# print an FST
%.txt: %.fst
	@fstprint --isymbols=$(SYMS) --osymbols=$(SYMS) $< > $@

################################################################################
# Miscellaneous                                                                #
################################################################################

clean:
	@rm -f fst/mutation.* fst/indel.* fst/marg_pos.* work/in_tape/* work/out_tape/* work/path/*
