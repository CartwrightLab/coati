# Construct an MG94 FST

library(stringr)
library(seqinr)
library(Matrix)

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

# construct codons and figure out which ones are synonymous
codons = cbind(rep(nucs,each=16),
               rep(nucs,times=4,each=4),
               rep(nucs,16))

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

# build the FST
r=1
for(i in 1:64) {
    for(j in 1:64) {
        arc(0,   r,   codons[i,1],codons[j,1],P[i,j])
        arc(r,   r+1, codons[i,2],codons[j,2])
        arc(r+1, 0,   codons[i,3],codons[j,3])
        r = r+2
    }
}

arc(0)
