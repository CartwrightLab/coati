# read the path through an FST and convert it to a FASTA alignment

library(data.table)
library(stringr)
library(seqinr)

# fix issues with seqinr
write.fasta = function (sequences, seqnames=names(sequences), file.out=stdout(), open = "w", nbchar = 60,
    as.string = FALSE)
{
    if(is.character(file.out)) {
    	outfile <- file(description = file.out, open = open)
    } else {
    	outfile = file.out
    }
    write.oneseq <- function(sequence, name, nbchar, as.string) {
        writeLines(paste(">", name, sep = ""), outfile)
        if (as.string)
            sequence <- s2c(sequence)
        l <- length(sequence)
        q <- floor(l/nbchar)
        r <- l - nbchar * q
        if (q > 0) {
            sapply(seq_len(q), function(x) writeLines(c2s(sequence[(nbchar *
                (x - 1) + 1):(nbchar * x)]), outfile))
        }
        if (r > 0) {
            writeLines(c2s(sequence[(nbchar * q + 1):l]), outfile)
        }
    }
    if (!is.list(sequences)) {
        write.oneseq(sequence = sequences, name = seqnames, nbchar = nbchar,
            as.string = as.string)
    }
    else {
        n.seq <- length(sequences)
        sapply(seq_len(n.seq), function(x) write.oneseq(sequence = as.character(sequences[[x]]),
            name = seqnames[x], nbchar = nbchar, as.string = as.string))
    }
    if(is.character(file.out)) {
	    close(outfile)
    }
}

# main function
fasta_main = function(input) {
	if(input == "-") {
		data = fread("cat /dev/stdin",fill=TRUE)
	} else {
		data = fread(file=input,fill=TRUE)
	}
	data[V3 == "<eps>", V3 := "-"]
	data[V4 == "<eps>", V4 := "-"]

	ret = list("seq_1"=data$V3, "seq_2"=data$V4)

	write.fasta(ret)
}

if(!interactive()) {
	ARGS = commandArgs(trailing=TRUE)
	fasta_main(input=ARGS[1])
}
