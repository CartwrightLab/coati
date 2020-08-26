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

#include <doctest/doctest.h>
#include <coati/align.hpp>

/* Alignment using dynamic programming implementation of marginal COATi model */
int mcoati(string fasta, vector<string> seq_names, vector<string> sequences,
	bool score, string weight_f, string output) {

	vector<string> alignment;
	float weight;
	ofstream out_w;

	if(score) {
		cout << alignment_score(sequences) << endl;
		return EXIT_SUCCESS;
	}

	alignment = mg94_marginal(sequences, weight);

	if(!weight_f.empty()) {
		// append weight and fasta file name to file
		out_w.open(weight_f, ios::app | ios::out);
		out_w << fasta << ",m-coati," << weight << endl;
		out_w.close();
	}

	// write alignment
	if(boost::filesystem::extension(output) == ".fasta") {
		return write_fasta(alignment, output, seq_names);
	} else {
		return write_phylip(alignment, output, seq_names);
	}
}

/* Alignment using FST library*/
int fst_alignment(string model, vector<VectorFst<StdArc>>& fsts, vector<string> seq_names,
	string fasta, string weight_f, string output, vector<string> sequences) {

	VectorFst<StdArc> mut_fst;

	if(model.compare("coati") == 0 ) {
		mg94(mut_fst);
	} else if(model.compare("dna")) {
		dna(mut_fst);
	} else if(model.compare("ecm")) {
		ecm(mut_fst);
	} else if(model.compare("m-ecm")) {
		ecm_marginal(mut_fst);
	} else {
		cout << "Mutation model unknown. Exiting!" << endl;
		return EXIT_FAILURE;
	}

	// get indel FST
	VectorFst<StdArc> indel_fst;
	indel(indel_fst, model);

	// sort mutation and indel FSTs
	VectorFst<StdArc> mutation_sort, indel_sort;
	mutation_sort = ArcSortFst<StdArc, OLabelCompare<StdArc>>(mut_fst, OLabelCompare<StdArc>());
	indel_sort = ArcSortFst<StdArc, ILabelCompare<StdArc>>(indel_fst, ILabelCompare<StdArc>());

	// compose mutation and indel FSTs
	ComposeFst<StdArc> coati_comp = ComposeFst<StdArc>(mutation_sort, indel_sort);

	// optimize coati FST
	VectorFst<StdArc> coati_fst;
	coati_fst = optimize(VectorFst<StdArc>(coati_comp));

	// find alignment graph
	// 1. compose in_tape and coati FSTs
	ComposeFst<StdArc> aln_inter = ComposeFst<StdArc>(fsts[0], coati_fst);
	// 2. sort intermediate composition
	VectorFst<StdArc> aln_inter_sort;
	aln_inter_sort = ArcSortFst<StdArc,OLabelCompare<StdArc>>(aln_inter, OLabelCompare<StdArc>());
	// 3. compose intermediate and out_tape FSTs
	VectorFst<StdArc> graph_fst;
	Compose(aln_inter_sort,fsts[1], &graph_fst);

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
		vector<StdArc::Weight> distance;
		ofstream out_w;

		ShortestDistance(aln_path, &distance);
		// append weight and fasta file name info in file
		out_w.open(weight_f, ios::app | ios::out);
		out_w << fasta << "," << model << "," << distance[0] << endl;
		out_w.close();
	}

	// topsort path FST
	TopSort(&aln_path);

	// write alignment
	if(boost::filesystem::extension(output) == ".fasta") {
		return write_fasta(aln_path, output, seq_names);
	} else {
		return write_phylip(aln_path, output, seq_names);
	}
}
