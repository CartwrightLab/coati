default: all

.PHONY: default all

all: fst/coati.fst

################################################################################
# Compiling                                                                    #
################################################################################

CXX = g++ -iquote src/include/coati/
CXXFLAGS_DEBUG = -Wall -g -lfst -O0 #-v
CXXFLAGS_EXEC = -lfst -Ofast
BOOSTFLAGS = -lboost_program_options

# -Wall turs on compiler warning flags
# -g build executable with debugging symbols
# -lfst link to shared FST library
# -O0 no optimization, faster compilation time
# -OFast high level of optimization, slower compile-time
# -iquote add directory to be search for header files when using #include "..."
# -v verbose

# compile
%.exe: src/lib/optimize.cc src/lib/matrix.cc src/lib/eigen.cc src/lib/mut_models.cc src/lib/utils.cc %.cc
	@$(CXX) $^ $(CXXFLAGS_DEBUG) $(BOOSTFLAGS) -o $@

%.fast: src/lib/optimize.cc src/lib/mut_models.cc src/lib/utils.cc %.cc
	@$(CXX) $^ $(CXXFLAGS_EXEC) $(BOOSTFLAGS) -o $@

# assembly code
%.s: %.cc
	@$(CXX) $(CXXFLAGS) -S -o $@ $<

################################################################################
# Construct raw FSTs from Rscripts. TODO: implement in c++                     #
# Construct fst/coati.fst                                                      #
# Original code: Reed A. Cartwright https://github.com/reedacartwright/toycoati#
################################################################################

RSCRIPT = Rscript --vanilla
# SYMS = fst/nuc_syms.txt
SYMS = fst/dna_syms.txt

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

# Construct an FST for a codon-marginalized MG94 codon model
fst/marg_pos.fst: scripts/marg_mutation.R $(SYMS)
	@$(RSCRIPT) $< 1 | fstcompile --arc_type=standard --isymbols=$(SYMS) \
		--osymbols=$(SYMS) - | fstrmepsilon | fstarcsort --sort_type=ilabel > $@

# Construct an FST for a nucleotide-marginalized MG94 codon model
fst/dna_marg.fst: scripts/marg_mutation.R $(SYMS)
	@$(RSCRIPT) $< 2 | fstcompile --arc_type=standard --isymbols=$(SYMS) \
		--osymbols=$(SYMS) - | fstrmepsilon | fstarcsort --sort_type=ilabel > $@

################################################################################
# Empirical Codon Model														   #
################################################################################
fst/ecm.fst: scripts/ecm.R $(SYMS)
	@$(RSCRIPT) $< | fstcompile --arc_type=standard --isymbols=$(SYMS) \
		--osymbols=$(SYMS) - | fstrmepsilon | fstarcsort --sort_type=ilabel > $@

################################################################################
# I/O acceptors & convert shortest path into an alignment                      #
# Original code: Reed A. Cartwright https://github.com/reedacartwright/toycoati#
################################################################################

# Extract input acceptor
work/in_tape/%.fst: fasta/% scripts/acceptor.R
	@$(RSCRIPT) scripts/acceptor.R $< 1 \
		| fstcompile --isymbols=$(SYMS) --osymbols=$(SYMS) - \
		| fstarcsort --sort_type=olabel > $@

# Extract output acceptor
work/out_tape/%.fst: fasta/% scripts/acceptor.R
	@$(RSCRIPT) scripts/acceptor.R $< 2 \
		| fstcompile --isymbols=$(SYMS) --osymbols=$(SYMS) - \
		| fstarcsort --sort_type=ilabel > $@

# Find alignment graph & its shortest path
#work/path/%.fst: src/main.exe work/in_tape/%.fst work/out_tape/%.fst
work/path/%.fst: build/coati work/in_tape/%.fst work/out_tape/%.fst
	@echo "Aligning "$*
	#@./src/main.fast -f $* -m dna -w aln/weights.csv # toycoati # marginalized
	@./build/coati -f $* -m ecm #-w aln/weights.csv # toycoati # marginalized

# Convert shortest path into an aligment
aln/%: work/path/%.fst scripts/fasta.R fst/indel.fst fst/marg_pos.fst fst/dna_marg.fst fst/mutation.fst
	@fstprint --isymbols=$(SYMS) --osymbols=$(SYMS) $< > work/path/$*.txt
	@cat work/path/$*.txt | $(RSCRIPT) scripts/fasta.R - > $@

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
# Testing                                                                      #
################################################################################

.PHONY: test

FASTA := $(shell ls ../toycoati/aln/)
ALN_FASTA = $(addprefix aln/,$(FASTA))
DIFF_FASTA = $(addsuffix .diff,$(FASTA))

# make sure fasta files in ../toycoati/aln/ are also in coati/fasta/
test: $(ALN_FASTA) $(DIFF_FASTA)

%.diff: aln/%
	diff aln/$* ../toycoati/aln/$*

################################################################################
# Miscellaneous                                                                #
################################################################################

clean:
	@rm -f fst/mutation.* fst/indel.* fst/marg_pos.* work/in_tape/* work/out_tape/* work/path/*
