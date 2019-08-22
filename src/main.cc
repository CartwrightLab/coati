
#include <fst/fstlib.h>
#include "lib/optimize.h"
#include <iostream>

using namespace fst;

int main(int argc, char *argv[]) {

	if(argc < 2){
		std::cout << "Name of fasta file to align needed." << std::endl;
		return -1;
	}

	std::string fasta = argv[1];

	// tropical semiring & read mutation and indel raw FSTs
	const VectorFst<StdArc> *mutation_raw = VectorFst<StdArc>::Read("fst/mutation.fst");
	const VectorFst<StdArc> *indel_raw = VectorFst<StdArc>::Read("fst/indel.fst");

	// optimize mutation and indel raw FSTs
	VectorFst<StdArc> mutation_fst, indel_fst;
	mutation_fst = optimize(*mutation_raw);
	indel_fst = optimize(*indel_raw);

	// sort mutation and indel FSTs
	VectorFst<StdArc> mutation_sort, indel_sort;
	mutation_sort = ArcSortFst<StdArc, OLabelCompare<StdArc>>(mutation_fst, OLabelCompare<StdArc>());
	indel_sort = ArcSortFst<StdArc, ILabelCompare<StdArc>>(indel_fst, ILabelCompare<StdArc>());

	// compose mutation and indel FSTs
	ComposeFst<StdArc> coati_comp = ComposeFst<StdArc>(mutation_sort, indel_sort);

	// optimize coati FST
	VectorFst<StdArc> coati_fst;
	coati_fst = optimize(VectorFst<StdArc>(coati_comp));

	// write coati FST
	// coati_fst.Write("fst/coati.fst");

	// read input and output FSAs
	const VectorFst<StdArc> *in_tape = VectorFst<StdArc>::Read("work/in_tape/"+fasta+".fst");
	const VectorFst<StdArc> *out_tape = VectorFst<StdArc>::Read("work/out_tape/"+fasta+".fst");

	std::cout << "Composing alignment graph..." << std::endl;
	// find alignment graph
	// 1. compose in_tape and coati FSTs
	ComposeFst<StdArc> aln_inter = ComposeFst<StdArc>(*in_tape, coati_fst);
	// 2. sort intermediate composition
	VectorFst<StdArc> aln_inter_sort;
	aln_inter_sort = ArcSortFst<StdArc,OLabelCompare<StdArc>>(aln_inter, OLabelCompare<StdArc>());
	// 3. compose intermediate and out_tape FSTs
	VectorFst<StdArc> graph_fst;
	Compose(aln_inter_sort,*out_tape, &graph_fst);

	//  case 1: ComposeFst is time: O(v1 v2 d1 (log d2 + m2)), space O(v1 v2)
	//			ShortestPath with ComposeFst (PDT) is time: O((V+E)^3), space: O((V+E)^3)
	//          Then convert to VectorFst (FST) is time, space: O(e^(O(V+E)))
	//  case 2: ComposeFst is time: O(v1 v2 d1 (log d2 + m2)), space: O(v1 v2)
	//			Convert ComposeFst (PDT) to VectorFst(FST) is time, space: O(e^(O(V+E)))
	//			Then ShortestPath with VectorFst is O(V log(V) + E)
	//  case 3: Compose is time: O(V1 V2 D1 (log D2 + M2)), space O(V1 V2 D1 M2)
	//			Then ShortestPath with VectorFst is O(V log(V) + E)

	std::cout << "Finding shortest path through graph..." << std::endl;
	// find shortest path through graph
	VectorFst<StdArc> aln_path;
	ShortestPath(graph_fst, &aln_path);

	std::cout << "Top sort..." << std::endl;
	// topsort path FST
	TopSort(&aln_path);

	// write path FST
	aln_path.Write("work/path/"+fasta+".fst");

    return 0;
}
