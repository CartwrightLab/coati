/*
# Copyright (c) 2020 Juan J. Garcia Mesa <juanjosegarciamesa@gmail.com>
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

/* Table for converting a nucleotide character to 2-bit encoding and for
	looking up coding amino acid based on codon index (AAA:0, AAC:1, ..., TTT:63) */
const uint8_t nt4_table[256] = {
	75, 78, 75, 78,  84, 84, 84, 84,  82, 83, 82, 83,  73, 73, 77, 73,
	81, 72, 81, 72,  80, 80, 80, 80,  82, 82, 82, 82,  76, 76, 76, 76,
	69, 68, 69, 68,  65, 65, 65, 65,  71, 71, 71, 71,  86, 86, 86, 86,
	42, 89, 42, 89,  83, 83, 83, 83,  42, 67, 87, 67,  76, 70, 76, 70,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

int read_fasta(string file, vector<string>& seq_names,
	vector<VectorFst<StdArc>>& fsts, vector<string>& sequences);
void add_arc(VectorFst<StdArc> &fst, int src, int dest, int ilabel=0,\
	int olabel=0, float weight=1.0);
VectorFst<StdArc> optimize(VectorFst<StdArc> fst);
int write_fasta(vector<string> alignment, string output, vector<string> seq_names);
int write_fasta(VectorFst<StdArc>& aln, string output, vector<string> seq_names);
int write_phylip(vector<string> alignment, string output, vector<string> seq_names);
int write_phylip(VectorFst<StdArc>& aln, string output, vector<string> seq_names);
bool acceptor(std::string content, VectorFst<StdArc> &accept);
int cod_distance(uint8_t cod1, uint8_t cod2);

#endif
