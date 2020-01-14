#include <fst/fstlib.h>
#include <iostream>
#include <boost/program_options.hpp>
#include <fstream>
#include "mut_models.h"
#include "utils.h"

using namespace fst;
namespace po = boost::program_options;

//TODO: incorporate unique pointers
int main(int argc, char *argv[]) {

	VectorFst<StdArc> mutation_fst;
	std::string fasta;
	std::string mut_model;
	std::string weight_f;

	try {
		po::options_description desc("Allowed options");
		desc.add_options()
			("help,h","arg1 fasta to align, arg2 model")
			("fasta,f",po::value<std::string>(&fasta)->required(), "name of fasta file")
			("model,m",po::value<std::string>(&mut_model)->required(), "substitution model")
			("weight,w",po::value<std::string>(&weight_f), "weight storing file")
		;

		po::variables_map varm;
		po::store(po::parse_command_line(argc,argv,desc),varm);

		if(varm.count("help") || varm.count("-h")) {
			std::cout << desc << std::endl;
			return 0;
		}

		po::notify(varm);

	} catch (po::error& e) {
		std::cerr << e.what() << ". Exiting!" << std::endl;
		return 4;
	}

	if(fasta.empty()) {
		std::cout << "fasta file required. Exiting!" << std::endl;
		return 1;
	}

	if(mut_model.empty()) {
	 	std::cout << "mutation model required. Exiting!" << std::endl;
	 	return 2;
	} else if(mut_model.compare("toycoati") == 0) {
		toycoati(mutation_fst);
	} else if(mut_model.compare("marginalized") == 0) {
		marg_mut(mutation_fst);
	} else if(mut_model.compare("dna") == 0) {
		dna_mut(mutation_fst);
	} else if(mut_model.compare("ecm") == 0) {
		ecm(mutation_fst);
	} else {
		std::cout << "Mutation model specified is unknown. Exiting!" << std::endl;
		return 3;
	}


	// TODO: think about deleting fst (delete fst;)

	// tropical semiring & read mutation and indel raw FSTs
	const VectorFst<StdArc> *indel_raw = VectorFst<StdArc>::Read("fst/indel.fst");

	// optimize mutation and indel raw FSTs
	VectorFst<StdArc> indel_fst;
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

	// read input and output FSAs
	const VectorFst<StdArc> *in_tape = VectorFst<StdArc>::Read("work/in_tape/"+fasta+".fst");
	const VectorFst<StdArc> *out_tape = VectorFst<StdArc>::Read("work/out_tape/"+fasta+".fst");

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

	// find shortest path through graph
	VectorFst<StdArc> aln_path;
	ShortestPath(graph_fst, &aln_path);

	// shortestdistance = weight of shortestpath
	if(!weight_f.empty()) {
		std::vector<StdArc::Weight> distance;
		std::ofstream out_w;

		ShortestDistance(aln_path, &distance);
		// append weight and fasta file name info in file
		out_w.open(weight_f, std::ios::app | std::ios::out);
		out_w << fasta << "," << mut_model << "," << distance[0] << std::endl;
		out_w.close();
	}

	// topsort path FST
	TopSort(&aln_path);

	// TODO: do aln_path -> fasta.fasta in C++
	// write path FST
	aln_path.Write("work/path/"+fasta+".fst");

    return 0;
}
