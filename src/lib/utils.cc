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

#include <coati/utils.hpp>

#define A nuc{'A',1,'R'}
#define C nuc{'C',2,'Y'}
#define G nuc{'G',3,'R'}
#define T nuc{'T',4,'Y'}
#define U nuc{'U',4,'Y'}
#define N nuc{'N',5,'-'}

#define PRINT_SIZE 100

/* comparison operator for nuc structure*/
bool operator==(nuc n1, nuc n2) {
	return(n1.nt == n2.nt);
}

/* comparison operator for cod structure*/
bool operator==(cod c1, cod c2) {
	return(c1.nt == c2.nt);
}

/* table of nuc[cleotides]*/
nuc nuc_table[6] = { A,C,G,T,U,N };

/* table of cod[ons]*/
cod cod_table[64] = {
	{{A,A,A},'K',1},{{A,A,C},'N',2},{{A,A,G},'K',3},{{A,A,T},'N',3},\
	{{A,C,A},'T',5},{{A,C,C},'T',6},{{A,C,G},'T',7},{{A,C,T},'T',8},\
	{{A,G,A},'R',9},{{A,G,C},'S',10},{{A,G,G},'R',11},{{A,G,T},'S',12},\
	{{A,T,A},'I',13},{{A,T,C},'I',14},{{A,T,G},'M',15},{{A,T,T},'I',16},\
	{{C,A,A},'Q',17},{{C,A,C},'H',18},{{C,A,G},'Q',19},{{C,A,T},'H',20},\
	{{C,C,A},'P',21},{{C,C,C},'P',22},{{C,C,G},'P',23},{{C,C,T},'P',24},\
	{{C,G,A},'R',25},{{C,G,C},'R',26},{{C,G,G},'R',27},{{C,G,T},'R',28},\
	{{C,T,A},'L',29},{{C,T,C},'L',30},{{C,T,G},'L',31},{{C,T,T},'L',32},\
	{{G,A,A},'E',33},{{G,A,C},'D',34},{{G,A,G},'E',35},{{G,A,T},'D',36},\
	{{G,C,A},'A',37},{{G,C,C},'A',38},{{G,C,G},'A',39},{{G,C,T},'A',40},\
	{{G,G,A},'G',41},{{G,G,C},'G',42},{{G,G,G},'G',43},{{G,G,T},'G',44},\
	{{G,T,A},'V',45},{{G,T,C},'V',46},{{G,T,G},'V',47},{{G,T,T},'V',48},\
	{{T,A,A},'*',49},{{T,A,C},'Y',50},{{T,A,G},'*',51},{{T,A,T},'Y',52},\
	{{T,C,A},'S',53},{{T,C,C},'S',54},{{T,C,G},'S',55},{{T,C,T},'S',56},\
	{{T,G,A},'*',57},{{T,G,C},'C',58},{{T,G,G},'W',59},{{T,G,T},'C',60},\
	{{T,T,A},'L',61},{{T,T,C},'F',62},{{T,T,G},'L',63},{{T,T,T},'F',64} };

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
int read_fasta(string file, vector<string>& seq_names, vector<VectorFst<StdArc>>& fsts) {
	ifstream input(file);
	if(!input.good()) {
		cerr << "Error opening '" << file << "'." << endl;
		exit(EXIT_FAILURE);
	}

	string line, name, content;
	while(getline(input, line).good() ) {
		if(line[0] == ';') continue;
		if(line.empty() || line[0] == '>') { // Identifier marker
			if(!name.empty()) {
				VectorFst<StdArc> accept;	// create FSA with sequence
				if(!acceptor(content, accept)) {
					cerr << "Creating acceptor from " << file << " failed. Exiting!" << endl;
					exit(EXIT_FAILURE);
				}
				fsts.push_back(accept);		// Add FSA
				name.clear();
			}
			if(!line.empty()) {		// Add name of sequence
				name = line.substr(1);
				seq_names.push_back(name);
			}
			content.clear();
		} else if(!name.empty()) {
			// might have to modify to delete spaces for simulation
			if(line.find(' ') == string::npos) {
				content += line;
			}
		}
	}
	if(!name.empty()) { // Add last sequence FSA if needed
		VectorFst<StdArc> accept;
		if(!acceptor(content, accept)) {
			cerr << "Creating acceptor from " << file << " failed. Exiting!" << endl;
			exit(EXIT_FAILURE);
		}
		fsts.push_back(accept);
	}

	return 0;
}

/* Write shortest path (alignment) in Fasta format */
void write_fasta(VectorFst<StdArc>& aln, string output, vector<string> seq_names) {
	SymbolTable *symbols = SymbolTable::ReadText("fst/dna_syms.txt");

	ofstream outfile;
	outfile.open(output);
	if(!outfile) {
		cerr << "Opening output file failed.\n";
		exit(EXIT_FAILURE);
	}

	string seq1, seq2;
	StdArc::StateId istate = aln.Start();
	StateIterator<StdFst> siter(aln);	// FST state iterator
	for(int i=0; i<aln.NumStates()-1; siter.Next(),i++) {
		ArcIterator<StdFst> aiter(aln, siter.Value());	// State arc iterator
		seq1.append(symbols->Find(aiter.Value().ilabel));
		seq2.append(symbols->Find(aiter.Value().olabel));
	}
	// map all epsilons (<eps>) to gaps (-)
	while(seq1.find("<eps>") != string::npos) seq1.replace(seq1.find("<eps>"),5,"-");
	while(seq2.find("<eps>") != string::npos) seq2.replace(seq2.find("<eps>"),5,"-");
	// write aligned sequences to file
	outfile << ">" << seq_names[0] << endl << seq1 << endl;
	outfile << ">" << seq_names[1] << endl << seq2 << endl;
	outfile.close();

}

/* Write shortest path (alignment) in PHYLIP format */
void write_phylip(VectorFst<StdArc>& aln, string output, vector<string> seq_names) {
	SymbolTable *symbols = SymbolTable::ReadText("fst/dna_syms.txt");

	ofstream outfile;
	outfile.open(output);
	if(!outfile) {
		cerr << "Opening output file failed.\n";
		exit(EXIT_FAILURE);
	}

	string seq1, seq2;
	StdArc::StateId istate = aln.Start();
	StateIterator<StdFst> siter(aln);	// FST state iterator
	for(int i=0; i<aln.NumStates()-1; siter.Next(),i++) {
		ArcIterator<StdFst> aiter(aln, siter.Value());	// State arc iterator
		seq1.append(symbols->Find(aiter.Value().ilabel));
		seq2.append(symbols->Find(aiter.Value().olabel));
	}
	// map all epsilons (<eps>) to gaps (-)
	while(seq1.find("<eps>") != string::npos) seq1.replace(seq1.find("<eps>"),5,"-");
	while(seq2.find("<eps>") != string::npos) seq2.replace(seq2.find("<eps>"),5,"-");

	assert(seq1.length() == seq2.length());

	// write aligned sequences to file
	outfile << seq_names.size() << " " << seq1.length() << endl;
	int i = PRINT_SIZE-4-max(seq_names[0].length(),seq_names[1].length());
	outfile << seq_names[0] << "\t" << seq1.substr(0,i) << endl;
	outfile << seq_names[1] << "\t" << seq2.substr(0,i) << endl << endl;
	for(; i<seq1.length(); i+=PRINT_SIZE) {
		outfile << seq1.substr(i,PRINT_SIZE) << endl;
		outfile << seq2.substr(i,PRINT_SIZE) << endl;
		outfile << endl;
	}

	outfile.close();

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
