#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <fstream>
#include <fst/fstlib.h>
#include <fst/script/print.h>
#include <fst/script/fst-class.h>
#include <fst/script/script-impl.h>

using namespace fst;
using namespace std;

/* nucleotide structure*/
struct nuc {
	char nt[1];	// nucleotide base
	int sym;	// fst symbol for that base
	char group;	// puRine or pYrimidine
};

/* codon structure*/
struct cod {
	nuc nt[3];		// codon nucleotides
	char subset;	// one of 20 subset groups according to Kioso et al. 2007
	int sym;		// fst symbol for the codon
};

/* nuc and cod structure comparison operators*/
bool operator==(nuc n1, nuc n2);
bool operator==(cod c1, cod c2);

extern nuc nuc_table[6];
extern cod cod_table[64];

int read_fasta(string file, vector<string>& seq_names, vector<VectorFst<StdArc>>& fsts);
void add_arc(VectorFst<StdArc> &fst, int src, int dest, int ilabel=0,\
	int olabel=0, float weight=1.0);
VectorFst<StdArc> optimize(VectorFst<StdArc> fst);
void write_fasta(VectorFst<StdArc>& aln, string output, vector<string> seq_names);
bool acceptor(std::string content, VectorFst<StdArc> &accept);

#endif
