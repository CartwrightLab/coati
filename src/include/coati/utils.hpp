/*
# Copyright (c) 2020 Juan J. Garcia Mesa <jgarc111@asu.edu>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
*/

#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <fstream>
#include <fst/fstlib.h>

using namespace fst;
using namespace std;

/* nucleotide structure*/
struct nuc {
	char nt[1];		// nucleotide base
	int sym;		// fst symbol for that base
	char group;		// puRine or pYrimidine
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
void write_phylip(VectorFst<StdArc>& aln, string output, vector<string> seq_names);
bool acceptor(std::string content, VectorFst<StdArc> &accept);

#endif
