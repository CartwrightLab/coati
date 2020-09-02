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

#include <fst/fstlib.h>
#include <iostream>
#include <boost/program_options.hpp>
#include <fstream>
#include <boost/filesystem.hpp>
#include <coati/mut_models.hpp>
#include <coati/align.hpp>

using namespace fst;
using namespace std;

namespace po = boost::program_options;

int main(int argc, char *argv[]) {

	string fasta, mut_model, weight_f, output, rate;
	bool score = false;

	try {
		po::options_description desc("Allowed options");
		desc.add_options()
			("help,h","Display this message")
			("fasta,f",po::value<string>(&fasta)->required(), "fasta file path")
			("model,m",po::value<string>(&mut_model)->default_value("m-coati"),
				"substitution model: coati, m-coati (default), dna, ecm, m-ecm")
			("weight,w",po::value<string>(&weight_f), "Write alignment score to file")
			("output,o",po::value<string>(&output), "Alignment output file")
			("score,s", "Calculate alignment score using marginal COATi model")
			("rate,r", po::value<string>(&rate), "Substitution rate matrix (CSV)")
		;

		po::positional_options_description pos_p;
		pos_p.add("fasta",-1);
		po::variables_map varm;
		po::store(po::command_line_parser(argc,argv).options(desc).positional(pos_p).run(),varm);

		if(varm.count("help")) {
			cout << desc << endl;
			return EXIT_SUCCESS;
		}

		if(varm.count("score")) {
			 score = true;
		}

		po::notify(varm);

	} catch (po::error& e) {
		cerr << e.what() << ". Exiting!" << endl;
		return EXIT_FAILURE;
	}

	// read input fasta file sequences as FSA (acceptors)
	vector<string> seq_names, sequences;
	vector<VectorFst<StdArc>> fsts;
	Matrix64f P;

	if(!rate.empty()) {
		Matrix64f Q;
		double br_len;
		parse_matrix_csv(rate, Q, br_len);
		// P matrix
		Q = Q * br_len;
		P = Q.exp();
	} else {
		P.setZero();
	}

	if(read_fasta(fasta,seq_names, fsts, sequences) != 0) {
		cerr << "Error reading " << fasta << " file. Exiting!" << endl;
		return EXIT_FAILURE;
	} else if(seq_names.size() < 2 || seq_names.size() != fsts.size()) {
		cerr << "At least two sequences required. Exiting!" << endl;
		return EXIT_FAILURE;
	}

	if(output.empty()) { // if no output is specified save in current dir in PHYLIP format
		output = boost::filesystem::path(fasta).stem().string()+".phy";
	} else if(boost::filesystem::extension(output) != ".phy" &&
		boost::filesystem::extension(output) != ".fasta") {
		cout << "Format for output file is not valid. Exiting!" << endl;
		return EXIT_FAILURE;
	}

	if(mut_model.compare("m-coati") == 0) {
		return mcoati(fasta, seq_names, sequences, score, weight_f, output, P);
	} else {
		return fst_alignment(mut_model, fsts, seq_names, fasta, weight_f, output, sequences);
	}
}
