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

#include <fst/fstlib.h>

#include <boost/program_options.hpp>
#include <coati/align.hpp>

namespace po = boost::program_options;

using namespace std;
using namespace fst;

int main(int argc, char* argv[]) {
    string rate;
    input_t in_data;

    try {
        po::options_description desc("Allowed options");
        desc.add_options()("help,h", "Display this message")(
            "fasta,f", po::value<std::filesystem::path>(&in_data.fasta_file.path)->required(),
            "fasta file path")(
            "model,m",
            po::value<string>(&in_data.mut_model)->default_value("m-coati"),
            "substitution model: coati, m-coati (default), dna, ecm, m-ecm")(
            "weight,w", po::value<string>(&in_data.weight_file),
            "Write alignment score to file")(
            "output,o", po::value<string>(&in_data.out_file),
            "Alignment output file")(
            "score,s",
            "Calculate alignment score using m-coati or m-ecm models")(
            "rate,r", po::value<string>(&in_data.rate),
            "Substitution rate matrix (CSV)")(
            "evo-time,t",
            po::value<double>(&in_data.br_len)->default_value(0.0133, "0.0133"),
            "Evolutionary time or branch length");

        po::positional_options_description pos_p;
        pos_p.add("fasta", -1);
        po::variables_map varm;
        po::store(po::command_line_parser(argc, argv)
                      .options(desc)
                      .positional(pos_p)
                      .run(),
                  varm);

        if(varm.count("help") || argc < 2) {
            cout << "Usage:	coati alignpair file.fasta [options]" << endl
                 << endl;
            cout << desc << endl;
            return EXIT_SUCCESS;
        }

        if(varm.count("score")) {
            in_data.score = true;
        }

        po::notify(varm);

    } catch(po::error& e) {
        cerr << e.what() << ". Exiting!" << endl;
        return EXIT_FAILURE;
    }

    // read input fasta file sequences as FSA (acceptors)
    vector<VectorFstStdArc> fsts;
    Matrix64f P;

    if(read_fasta(in_data.fasta_file, fsts) != 0) {
        cerr << "Error reading " << in_data.fasta_file.path << " file. Exiting!"
             << endl;
        return EXIT_FAILURE;
    } else if(in_data.fasta_file.seq_names.size() != 2 ||
              in_data.fasta_file.seq_names.size() != fsts.size()) {
        cerr << "Exactly two sequences required. Exiting!" << endl;
        return EXIT_FAILURE;
    }

    if(in_data.out_file.empty()) {  // if no output is specified save in current
                                    // dir in PHYLIP format
        in_data.out_file =
            std::filesystem::path(in_data.fasta_file.path).stem().string() +
            ".phy";
    } else if(std::filesystem::path(in_data.out_file).extension() != ".phy" &&
              std::filesystem::path(in_data.out_file).extension() != ".fasta") {
        cout << "Format for output file is not valid. Exiting!" << endl;
        return EXIT_FAILURE;
    }

    if(!rate.empty()) {
        in_data.mut_model = "user_marg_model";

        Matrix64f Q;
        double br_len;
        parse_matrix_csv(rate, Q, br_len);
        // P matrix
        Q = Q * br_len;
        P = Q.exp();

        return mcoati(in_data, P);
    } else if((in_data.mut_model.compare("m-coati") == 0) ||
              in_data.mut_model.compare("no_frameshifts") == 0) {
        mg94_p(P, in_data.br_len);
        return mcoati(in_data, P);
    } else if(in_data.mut_model.compare("m-ecm") == 0) {
        ecm_p(P, in_data.br_len);
        return mcoati(in_data, P);
    } else {
        return fst_alignment(in_data, fsts);
    }
}
