/*
# Copyright (c) 2020-2021 Juan J. Garcia Mesa <juanjosegarciamesa@gmail.com>
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
#include <coati/dna_syms.hpp>

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
	SymbolTable syms;
	fill_symbol_table(syms);

	EncodeMapper<StdArc> encoder(kEncodeLabels, ENCODE);
	encoder.SetInputSymbols(&syms);
	encoder.SetOutputSymbols(&syms);
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
int read_fasta(fasta_t& fasta_file, vector<VectorFst<StdArc>>& fsts) {
	ifstream input(fasta_file.path);
	if(!input.good()) {
		cerr << "Error opening '" << fasta_file.path << "'." << endl;
		return EXIT_FAILURE;
	}

	string line, name, content;
	while(getline(input, line).good() ) {
		if(line[0] == ';') continue;
		if(line[0] == '>') { // Identifier marker
			if(!name.empty()) {
				VectorFst<StdArc> accept;	// create FSA with sequence
				if(!acceptor(content, accept)) {
					cerr << "Creating acceptor from " << fasta_file.path <<
						" failed. Exiting!" << endl;
					return EXIT_FAILURE;
				}
				fsts.push_back(accept);		// Add FSA
				fasta_file.seq_data.push_back(content);
				name.clear();
			}
			// Add name of sequence
			name = line.substr(1);
			fasta_file.seq_names.push_back(name);
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
			cerr << "Creating acceptor from " << fasta_file.path << " failed. Exiting!" << endl;
			return EXIT_FAILURE;
		}
		fsts.push_back(accept);
		fasta_file.seq_data.push_back(content);
	}

	return 0;
}

/* Read fasta format file */
int read_fasta(fasta_t& fasta_file) {
	ifstream input(fasta_file.path);
	if(!input.good()) {
		cerr << "Error opening '" << fasta_file.path << "'." << endl;
		return EXIT_FAILURE;
	}

	string line, name, content;
	while(getline(input, line).good() ) {
		if(line[0] == ';') continue;
		if(line[0] == '>') { // Identifier marker
			if(!name.empty()) {
				fasta_file.seq_data.push_back(content);
				name.clear();
			}
			// Add name of sequence
			name = line.substr(1);
			fasta_file.seq_names.push_back(name);
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
		fasta_file.seq_data.push_back(content);
	}

	return 0;
}

/* Write alignment in Fasta format */
int write_fasta(fasta_t& fasta_file) {
	ofstream outfile;
	outfile.open(fasta_file.path);
	if(!outfile) {
		cerr << "Opening output file failed.\n";
		return EXIT_FAILURE;
	}

	for(int i=0; i < fasta_file.seq_names.size(); i++) {
		outfile << ">" << fasta_file.seq_names[i] << endl << fasta_file.seq_data[i] << endl;
	}
	outfile.close();

	return EXIT_SUCCESS;
}

/* Write shortest path (alignment) in Fasta format */
int write_fasta(VectorFst<StdArc>& aln, fasta_t& fasta_file) {
	SymbolTable symbols;
	fill_symbol_table(symbols);

	string seq1, seq2;
	StateIterator<StdFst> siter(aln);	// FST state iterator
	for(int i=0; i<aln.NumStates()-1; siter.Next(),i++) {
		ArcIterator<StdFst> aiter(aln, siter.Value());	// State arc iterator
		seq1.append(symbols.Find(aiter.Value().ilabel));
		seq2.append(symbols.Find(aiter.Value().olabel));
	}

	fasta_file.seq_data.push_back(seq1);
	fasta_file.seq_data.push_back(seq2);

	// map all epsilons (<eps>) to gaps (-)
	while(fasta_file.seq_data[0].find("<eps>") != string::npos) fasta_file.seq_data[0].replace(fasta_file.seq_data[0].find("<eps>"),5,"-");
	while(fasta_file.seq_data[1].find("<eps>") != string::npos) fasta_file.seq_data[1].replace(fasta_file.seq_data[1].find("<eps>"),5,"-");

	// write aligned sequences to file
	write_fasta(fasta_file);
	return EXIT_SUCCESS;
}

/* Write alignment in PHYLIP format */
int write_phylip(fasta_t& fasta_file) {
	ofstream outfile;
	outfile.open(fasta_file.path);
	if(!outfile) {
		cerr << "Opening output file failed.\n";
		return EXIT_FAILURE;
	}

	// write aligned sequences to file
	outfile << fasta_file.seq_names.size() << " " << fasta_file.seq_data[0].length() << endl;
	int i = PRINT_SIZE-4-max(fasta_file.seq_names[0].length(),fasta_file.seq_names[1].length());
	for(int j=0; j < fasta_file.seq_names.size(); j++) {
		outfile << fasta_file.seq_names[j] << "\t" << fasta_file.seq_data[j].substr(0,i) << endl;
	}
	outfile << endl;

	for(; i<fasta_file.seq_data[0].length(); i+=PRINT_SIZE) {
		for(int j=0; j < fasta_file.seq_names.size(); j++) {
			outfile << fasta_file.seq_data[j].substr(i,PRINT_SIZE) << endl;
		}
		outfile << endl;
	}

	return EXIT_SUCCESS;
}

/* Write shortest path (alignment) in PHYLIP format */
int write_phylip(VectorFst<StdArc>& aln, fasta_t& fasta_file) {
	SymbolTable symbols;
	fill_symbol_table(symbols);

	string seq1, seq2;
	StateIterator<StdFst> siter(aln);	// FST state iterator
	for(int i=0; i<aln.NumStates()-1; siter.Next(),i++) {
		ArcIterator<StdFst> aiter(aln, siter.Value());	// State arc iterator
		seq1.append(symbols.Find(aiter.Value().ilabel));
		seq2.append(symbols.Find(aiter.Value().olabel));
	}

	fasta_file.seq_data.push_back(seq1);
	fasta_file.seq_data.push_back(seq2);

	// map all epsilons (<eps>) to gaps (-)
	while(fasta_file.seq_data[0].find("<eps>") != string::npos) fasta_file.seq_data[0].replace(fasta_file.seq_data[0].find("<eps>"),5,"-");
	while(fasta_file.seq_data[1].find("<eps>") != string::npos) fasta_file.seq_data[1].replace(fasta_file.seq_data[1].find("<eps>"),5,"-");

	write_phylip(fasta_file);
	return EXIT_SUCCESS;
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

/* Cast codon to position in codon list AAA->0, AAAC->1 ... TTT->63 */
int cod_int(string codon) {
	return ((uint8_t) nt4_table[codon[0]]<<4)+((uint8_t) nt4_table[codon[1]]<<2)
		+((uint8_t) nt4_table[codon[2]]);
}


/* Read substitution rate matrix from a CSV file */
int parse_matrix_csv(string file, Matrix64f& P, double& br_len) {
	ifstream input(file);
	if(!input.good()) {
		cerr << "Error opening '" << file << "'." << endl;
		return EXIT_FAILURE;
	}

	string line;
	if(input.good()) {
		// Read branch length
		getline(input,line);
		br_len = stod(line);
	}

	vector<string> vec;
	int count = 0;

	while(getline(input, line)) {
		boost::algorithm::split(vec,line,boost::is_any_of(","));
		P(cod_int(vec[0]), cod_int(vec[1])) = stod(vec[2]);
		count++;
	}

	if(count != 4096){
		cout << "Error reading substitution rate CSV file. Exiting!" << endl;
		return EXIT_FAILURE;
	}

	return 0;
}
