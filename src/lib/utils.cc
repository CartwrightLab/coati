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
#include <coati/utils.hpp>

#define PRINT_SIZE 100

/* Add arc to FST */
void add_arc(VectorFst<StdArc> &n2p, int src, int dest, int ilabel, int olabel, float weight) {
	if(weight == 1.0) {
		weight = 0.0;
	} else if(weight == 0.0) {
		weight = INT_MAX;
	} else {
		weight = -log(weight);
	}

	if(n2p.NumStates() <= dest) {
		n2p.AddState();
		n2p.AddArc(src, StdArc(ilabel, olabel, weight, dest));
	} else {
		n2p.AddArc(src, StdArc(ilabel, olabel, weight, dest));
	}
}

/* Optimize FST: remove epsilons, determinize, and minimize */
VectorFst<StdArc> optimize(VectorFst<StdArc> fst_raw) {
	// encode FST
	SymbolTable *syms = SymbolTable::ReadText("fst/dna_syms.txt");
	EncodeMapper<StdArc> encoder(kEncodeLabels, ENCODE);
	encoder.SetInputSymbols(syms);
	encoder.SetOutputSymbols(syms);
	EncodeFst<StdArc> fst_enc = EncodeFst<StdArc>(fst_raw, encoder);

	// reduce to more efficient form
	// 1. epsilon removal
	VectorFst<StdArc> fst_rmep;
	fst_rmep = RmEpsilonFst<StdArc>(fst_enc);

	// 2. determinize
	VectorFst<StdArc> fst_det;
	fst_det = DeterminizeFst<StdArc>(fst_rmep);

	// 3. minimize
	VectorFst<StdArc> fst_min;
	Minimize(&fst_det, &fst_min);

	// decode
	DecodeFst<StdArc> fst_dec = DecodeFst<StdArc>(fst_det, encoder);

	// epsilon removal after decoding
	VectorFst<StdArc> fst_opt;
	fst_opt = RmEpsilonFst<StdArc>(fst_dec);

	return fst_opt;
}

/* Read fasta format file */
int read_fasta(string file, vector<string>& seq_names,
		vector<VectorFst<StdArc>>& fsts, vector<string>& sequences) {
	ifstream input(file);
	if(!input.good()) {
		cerr << "Error opening '" << file << "'." << endl;
		exit(EXIT_FAILURE);
	}

	string line, name, content;
	while(getline(input, line).good() ) {
		if(line[0] == ';') continue;
		if(line[0] == '>') { // Identifier marker
			if(!name.empty()) {
				VectorFst<StdArc> accept;	// create FSA with sequence
				if(!acceptor(content, accept)) {
					cerr << "Creating acceptor from " << file << " failed. Exiting!" << endl;
					exit(EXIT_FAILURE);
				}
				fsts.push_back(accept);		// Add FSA
				sequences.push_back(content);
				name.clear();
			}
			// Add name of sequence
			name = line.substr(1);
			seq_names.push_back(name);
			content.clear();
		} else if(line.empty()) {
			continue;	// omit empty lines
		} else if(!name.empty()) {
			// Remove spaces
			line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());
			content += line;
		}
	}
	if(!name.empty()) { // Add last sequence FSA if needed
		VectorFst<StdArc> accept;
		if(!acceptor(content, accept)) {
			cerr << "Creating acceptor from " << file << " failed. Exiting!" << endl;
			exit(EXIT_FAILURE);
		}
		fsts.push_back(accept);
		sequences.push_back(content);
	}

	return 0;
}

/* Write alignment in Fasta format */
void write_fasta(vector<string> alignment, string output, vector<string> seq_names) {
	ofstream outfile;
	outfile.open(output);
	if(!outfile) {
		cerr << "Opening output file failed.\n";
		exit(EXIT_FAILURE);
	}

	outfile << ">" << seq_names[0] << endl << alignment[0] << endl;
	outfile << ">" << seq_names[1] << endl << alignment[1] << endl;
	outfile.close();

}

/* Write shortest path (alignment) in Fasta format */
void write_fasta(VectorFst<StdArc>& aln, string output, vector<string> seq_names) {
	SymbolTable *symbols = SymbolTable::ReadText("fst/dna_syms.txt");

	vector<string> alignment;// seq1, seq2;
	string seq1, seq2;
	StdArc::StateId istate = aln.Start();
	StateIterator<StdFst> siter(aln);	// FST state iterator
	for(int i=0; i<aln.NumStates()-1; siter.Next(),i++) {
		ArcIterator<StdFst> aiter(aln, siter.Value());	// State arc iterator
		seq1.append(symbols->Find(aiter.Value().ilabel));
		seq2.append(symbols->Find(aiter.Value().olabel));
	}

	alignment.push_back(seq1);
	alignment.push_back(seq2);

	// map all epsilons (<eps>) to gaps (-)
	while(alignment[0].find("<eps>") != string::npos) alignment[0].replace(alignment[0].find("<eps>"),5,"-");
	while(alignment[1].find("<eps>") != string::npos) alignment[1].replace(alignment[1].find("<eps>"),5,"-");

	// write aligned sequences to file
	write_fasta(alignment, output, seq_names);

}

/* Write alignment in PHYLIP format */
void write_phylip(vector<string> alignment, string output, vector<string> seq_names) {
	ofstream outfile;
	outfile.open(output);
	if(!outfile) {
		cerr << "Opening output file failed.\n";
		exit(EXIT_FAILURE);
	}

	// write aligned sequences to file
	outfile << seq_names.size() << " " << alignment[0].length() << endl;
	int i = PRINT_SIZE-4-max(seq_names[0].length(),seq_names[1].length());
	outfile << seq_names[0] << "\t" << alignment[0].substr(0,i) << endl;
	outfile << seq_names[1] << "\t" << alignment[1].substr(0,i) << endl << endl;
	for(; i<alignment[0].length(); i+=PRINT_SIZE) {
		outfile << alignment[0].substr(i,PRINT_SIZE) << endl;
		outfile << alignment[1].substr(i,PRINT_SIZE) << endl;
		outfile << endl;
	}


}

/* Write shortest path (alignment) in PHYLIP format */
void write_phylip(VectorFst<StdArc>& aln, string output, vector<string> seq_names) {
	SymbolTable *symbols = SymbolTable::ReadText("fst/dna_syms.txt");

	vector<string> alignment;
	string seq1, seq2;
	StdArc::StateId istate = aln.Start();
	StateIterator<StdFst> siter(aln);	// FST state iterator
	for(int i=0; i<aln.NumStates()-1; siter.Next(),i++) {
		ArcIterator<StdFst> aiter(aln, siter.Value());	// State arc iterator
		seq1.append(symbols->Find(aiter.Value().ilabel));
		seq2.append(symbols->Find(aiter.Value().olabel));
	}

	alignment.push_back(seq1);
	alignment.push_back(seq2);

	// map all epsilons (<eps>) to gaps (-)
	while(alignment[0].find("<eps>") != string::npos) alignment[0].replace(alignment[0].find("<eps>"),5,"-");
	while(alignment[1].find("<eps>") != string::npos) alignment[1].replace(alignment[1].find("<eps>"),5,"-");

	assert(alignment[0].length() == alignment[1].length());

	write_phylip(alignment, output, seq_names);

}

// TEST_CASE("[utils.cc] write_phylip") {
//
// }

/* Create FSAs (acceptors) from a fasta file*/
bool acceptor(string content, VectorFst<StdArc> &accept) {
	map<int, char> syms = \
		{{'-',0},{'A',1},{'C',2},{'G',3},{'T',4},{'U',4},{'N',5}};

	// Add initial state
	accept.AddState();
	accept.SetStart(0);

	for(int i=0; i<content.length(); i++) {
		add_arc(accept, i, i+1, syms[content.at(i)], syms[content.at(i)]);
	}

	// Add final state and run an FST sanity check (Verify)
	accept.SetFinal(accept.NumStates()-1,0.0);
	return Verify(accept);
}

/* Hamming distance between two codons */
int cod_distance(uint8_t cod1, uint8_t cod2) {
	int distance = 0;

	distance += (((cod1 & 48) >> 4) == ((cod2 & 48) >> 4) ? 0:1);
	distance += (((cod1 & 12) >> 2) == ((cod2 & 12) >> 2) ? 0:1);
	distance += ((cod1 & 3) == (cod2 & 3) ? 0:1);

	return distance;
}
