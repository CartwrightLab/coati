#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <fstream>
#include <fst/fstlib.h>

extern std::map<int, char> nuc_sym;
extern char nuc_table[6];

void add_arc(fst::VectorFst<fst::StdArc> &fst, int src, int dest, int ilabel=0,\
	int olabel=0, float weight=1.0);
fst::VectorFst<fst::StdArc> optimize(fst::VectorFst<fst::StdArc>);
void write_fasta(fst::VectorFst<fst::StdArc>& aln, std::string fasta, std::string outdir);

#endif
