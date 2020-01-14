#include <iostream>
#include <fstream>
#include <fst/fstlib.h>

using namespace fst;

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

// TODO: write as template for more fst and arc type variability
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
