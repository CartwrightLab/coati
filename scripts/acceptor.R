# convert a sequence in a fasta file into an acceptor
# Output format is described at http://www.openfst.org/twiki/bin/view/FST/FstQuickTour#CreatingFsts
 
library(stringr)
library(seqinr)

block_width = 1

acceptor_main = function(input, species) {
	seqs = read.fasta(input, set.attributes=FALSE)
	a = seqs[[species]] %>% toupper()
	cat(sprintf("%d\t%d\t%s\t%s\n",seq_along(a),seq_along(a)+1,a,a), sep="")
	cat(sprintf("%d\n", length(a)+1))

}

if(!interactive()) {
	ARGS = commandArgs(trailing=TRUE)
	acceptor_main(input=ARGS[1],species=ARGS[2])
}
