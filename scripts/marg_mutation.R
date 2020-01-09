
library(Matrix)
library(stringr)
library(seqinr)


get_codons = function() {
  nucs = c("A","C","G","T")
  codons = cbind(rep(nucs,each=16),
                 rep(nucs,times=4,each=4),
                 rep(nucs,16))
  return(codons)
}

# original code from mutation.R - toycoati - Reed Cartwright
MG94 = function(){
  # parameters
  # Yang (1994) Estimating the pattern of nucleotide substitution
  nucs = c("A","C","G","T")

  nuc_freqs = c(0.308,0.185,0.199,0.308)

  nuc_q = c(-0.818, 0.132, 0.586, 0.1,
            0.221, -1.349, 0.231, 0.897,
            0.909, 0.215, -1.322, 0.198,
            0.1, 0.537, 0.128, -0.765)
  nuc_q = matrix(nuc_q,4,4,byrow=T)

  omega = 0.2
  brlen = 0.0133


  # construct codons and figure out which ones are synonymous
  codons = get_codons()
  codonstrs = apply(codons,1,str_c,collapse="")

  syn = syncodons(codonstrs)
  names(syn) = toupper(names(syn))
  syn = lapply(syn,toupper)

  # MG94 model - doi:10.1534/genetics.108.092254
  Q = matrix(0,64,64)
  Pi = rep(0,64)
  # construct the transition matrix
  for(i in 1:64) {
    Pi[i] = prod(nuc_freqs[match(codons[i,],nucs)])

    for(j in 1:64) {
      if(i == j) {
        Q[i,j] = 0
      } else if(sum(codons[i,] != codons[j,]) > 1) {
        Q[i,j] = 0
      } else {
        if(codonstrs[j] %in% syn[[ codonstrs[i] ]]) {
          w = 1
        } else {
          w = omega
        }
        o = which(codons[i,] != codons[j,])
        x = which(nucs == codons[i,o])
        y = which(nucs == codons[j,o])

        Q[i,j] = w*nuc_q[x,y];
      }
    }
  }

  # normalize Q
  diag(Q) = -rowSums(Q)
  Q = Q / -sum(Pi*diag(Q))

  # construct transition matrix
  P = expm(Q*brlen)
  return(P)
}

# Function for creating a transition in an FST
# Output format is described at http://www.openfst.org/twiki/bin/view/FST/FstQuickTour#CreatingFsts
arc = function(src, dest=NULL, ilab="<eps>", olab="<eps>", weight=1.0) {
  if(weight == 1.0) {
    weight = 0.0
  } else {
    weight = -log(weight)
  }

  if(is.null(dest)) {
    cat(sprintf("%i\t%g\n", src, weight))
  } else {
    cat(sprintf("%i\t%i\t%s\t%s\t%g\n", src, dest, ilab, olab, weight))
  }
}

# # create FST to read 3 nucleotides and output a codon num (1-64)
# # TODO: codonstrs = apply(codons,1,str_c,collapse="")
# nuc2cod = function() {
# 	codons = get_codons()
# 	s = 1
# 	for(i in 1:nrow(get_codons())) {
# 		arc(0,		s,		codons[i,1],	"<eps>")
# 		arc(s,		s+1,	codons[i,2],	"<eps>")
# 		arc(s+1,	0,		codons[i,3],	paste(codons[i,],collapse="",sep = ""))
# 		s = s+2
# 	}
# 	arc(0)
#
# }

# # create FST to read 1 codon and output 3 (1-192) positions
# cod2pos = function() {
# 	codons = get_codons()
# 	s = 1
# 	for(i in 1:64) {
# 		arc(0,		s,		paste(codons[i,],collapse="",sep=""),	paste(paste(codons[i,],collapse="",sep = ""),1,sep = ""))
# 		arc(s,		s+1,	"<eps>",								paste(paste(codons[i,],collapse="",sep = ""),2,sep = ""))
# 		arc(s+1,	0,		"<eps>",								paste(paste(codons[i,],collapse="",sep = ""),3,sep = ""))
# 		s = s+2
# 	}
# 	arc(0)
# }

# nuc2pos = function() {
# 	nucs = c("A","C","G","T")
# 	codons = get_codons()
# 	codonstrs = apply(codons,1,str_c,collapse="")
# 	s = o = 1
# 	for(i in 1:4) {
# 		for(j in 1:4) {
# 			for(k in 1:4) {
# 				arc(0,		s,		nucs[i],	paste(codonstrs[o],1,collapse="",sep=""))
# 				arc(s,		s+1,	nucs[j],	paste(codonstrs[o],2,collapse="",sep=""))
# 				arc(s+1,	0,		nucs[k],	paste(codonstrs[o],3,collapse="",sep=""))
# 				s = s+2
# 				o = o+1
# 			}
# 		}
# 	}
# 	arc(0)
# }

marg_pos = function() {
	P = MG94()
	codons = get_codons()
	nucs = c("A","C","G","T")

	codons = get_codons()
  	codonstrs = apply(codons,1,str_c,collapse="")

	p_marg = array(0, dim = c(64,3,4))
	for(i in 1:64) {
		for(j in 1:3) {
			for(k in 1:4) {
				p_marg[i,j,k] = sum(P[i,codons[,j]==nucs[k]])
				arc(0,	0,	paste(codonstrs[i],j,sep=""),	nucs[k],	p_marg[i,j,k])
			}
		}
	}

	# for(i in 1:64) {
	# 	for(j in 1:3) {
	# 		for(k in 1:4) {
	# 			arc(0,	0,	paste(codonstrs[i],j,sep=""),	nucs[k],	p_marg[i,j,k])
	# 		}
	# 	}
	# }
	arc(0)
}

dna_marg = function() {
	# why redo code from function above?
	P = MG94()
	codons = get_codons()
	nucs = c("A","C","G","T")

	codons = get_codons()
  	codonstrs = apply(codons,1,str_c,collapse="")

	dna_m = matrix(0,4,4)

	for(i in 1:64) {  # for 1st row
		for(j in 1:3) {   # for position 1
			for(k in 1:4) {   # for from A
				for(l in 1:4) {   # for to A
					dna_m[k,l] = dna_m[k,l] + sum(P[i,codons[,j]==nucs[l]]*ifelse(codons[i,j]==nucs[k],1,0))
					# maybe I'm counting extra and only need to count each codon to codon once if both contain a "match" (N_i=nuc & M_i=nuc)
				}
			}
		}
	}

	dna_m = dna_m / rowSums(dna_m)

	for(i in 1:4) {
		for(j in 1:4) {
			arc(0,	0,	nucs[i],	nucs[j],	dna_m[i,j])
		}
	}
	arc(0)
	# return(dna_m)
}


if(!interactive()) {
    ARGS = commandArgs(trailing=TRUE)
	if(length(ARGS) < 1) {
		print(paste("At least one integer argument needed:"))
		print(paste("	1 (marg_mut), 2 (dna_marg)"))
	}
	else {
		# TODO: improve sloppy code!
		for(i in length(ARGS)) {
			if(ARGS[i] == 1) {marg_pos()}
			if(ARGS[i] == 2) {dna_marg()}
		}
	}
}
