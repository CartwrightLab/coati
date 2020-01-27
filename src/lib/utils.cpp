#include "utils.h"

using namespace fst;

#define A nuc{'A',1,'R'}
#define C nuc{'C',2,'Y'}
#define G nuc{'G',3,'R'}
#define T nuc{'T',4,'Y'}
#define U nuc{'U',4,'Y'}
#define N nuc{'N',5,'-'}

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
		weight = -std::log(weight);
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

// TODO: change to int type and update fst by reference to check exit status
int read_fasta(std::string file, std::string fasta) {
	VectorFst<StdArc> acceptor;

	std::ifstream input(file);
	if(!input.good()) {
		std::cerr << "Error opening '" << file << "'." << std::endl;
	}

	std::string line, name, content;
	while(std::getline(input, line).good() ) {
		if(!line.empty() && line[0] == ';') {
			; // do nothing if comment line or empty line
		} else if(line.empty()) {
			if(content.empty()) {
				return -1;
			} else {
				content.clear();
			}
		} else if(line[0] == '>') { // Identifier marker
			name = line.substr(1);
			content.clear();
		} else if(!name.empty()) {
			content += line;
		}
	}



	return 0;
}

// TODO: read and write proper sequence ID
/* Write shortest path (alignment) in Fasta format */
void write_fasta(VectorFst<StdArc>& aln, std::string fasta, std::string outdir) {
	std::map<int, char> syms = \
		{{0,'-'},{1,'A'},{2,'C'},{3,'G'},{4,'T'},{4,'U'},{5,'N'}};

	std::ofstream outfile;
	outfile.open(outdir+"/"+fasta);
	if(!outfile) {
		std::cout << "Opening output file failed.";
		exit(EXIT_FAILURE);
	}

	outfile << "Seq_1" << std::endl;
	int aln_len = aln.NumStates()-1;
	char seq_2[aln_len];

	StdArc::StateId istate = aln.Start();
	StateIterator<StdFst> siter(aln);
	for(int i=0; i<aln_len; siter.Next(),i++) {
		ArcIterator<StdFst> aiter(aln, siter.Value());
		outfile << syms[aiter.Value().ilabel];
		seq_2[i] = syms[aiter.Value().olabel];
	}
	outfile << std::endl << "Seq_2" << std::endl;
	outfile.write(seq_2,aln_len);
	outfile << std::endl;
	outfile.close();
}

// int acceptor(std::string file, VectorFst<StdArc> &acceptor) {
//
// 	return 0;
// }

/*
// Original code from EdwardsLab: https://edwards.sdsu.edu/research/fastq-to-fasta/
// convert a fastq file to fasta
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
using namespace std;

int main (int argc, char* argv[]) {

  if ( argc < 3) {
    cerr << "Usage: " << argv[0] << " <fastq file (use - to read from STDIN)> <fasta file>\n";
    return 1;
  }


  ifstream fastq;
  streambuf* orig_cin = 0;
  if (strcmp(argv[1], "-") != 0) {
    cout << "reading from " << argv[1] << '\n';
    fastq.open(argv[1]);
    if (!fastq) return 1;
    orig_cin = cin.rdbuf(fastq.rdbuf());
    cin.tie(0); // tied to cout by default
  }

  string line;
  ofstream fasta(argv[2]);
  if (!fasta) return 2;
  int c=0;

  while (getline(cin, line))
  {
    //cout << c << " : " << line << '\n';
    if ( c==0 ) {
      line.replace(0, 1, ">");
      fasta << line << '\n';
    }
    if ( c==1 ) {
      fasta << line << '\n';
    }
    c++;
    if ( c == 4) { c=0 ;}
  }

  if ( c != 0 ) {
	  cerr << "ERROR: There appears to be the wrong number of lines in your file!" << endl;
	return 1;
  }

  return 0;
}
*/
