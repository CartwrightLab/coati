#ifndef UTILS_H
#define UTILS_H

#include <fst/fstlib.h>

void add_arc(fst::VectorFst<fst::StdArc> &fst, int src, int dest, int ilabel=0,\
	int olabel=0, float weight=1.0);
// void add_arc(fst::VectorFst<fst::StdArc> &n2p, int src, int dest, int ilabel, int olabel, float weight);
// int read2Fasta(std::string file);

#endif
