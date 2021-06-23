/*
# Copyright (c) 2021 Juan J. Garcia Mesa <juanjosegarciamesa@gmail.com>
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

#include <boost/program_options.hpp>
#include <coati/align.hpp>

namespace po = boost::program_options;

int main(int argc, char* argv[]) {
    input_t in_data;

    try {
        po::options_description desc("Allowed options");
        desc.add_options()("help,h", "Display this message")(
            "fasta,f", po::value<string>(&in_data.fasta_file.path)->required(),
            "fasta file path")(
            "model,m",
            po::value<string>(&in_data.mut_model)->default_value("m-coati"),
            "substitution model: m-coati (default), m-ecm")(
            "weight,w", po::value<string>(&in_data.weight_file),
            "write alignment score to given file")(
            "output,o", po::value<string>(&in_data.out_file),
            "alignment output file")(
            "tree,p", po::value<string>(&in_data.tree)->required(),
            "newick phylogenetic tree")(
            "ref,r", po::value<string>(&in_data.ref)->required(),
            "reference sequence");

        po::positional_options_description pos_p;
        pos_p.add("fasta", 1);
        pos_p.add("tree", 1);
        po::variables_map varm;
        po::store(po::command_line_parser(argc, argv)
                      .options(desc)
                      .positional(pos_p)
                      .run(),
                  varm);

        if(varm.count("help") || argc < 3) {
            cout << "Usage: coati msa file.fasta tree.newick [options]" << endl
                 << endl;
            cout << desc << endl;
            return EXIT_SUCCESS;
        }

        po::notify(varm);

    } catch(po::error& e) {
        cerr << e.what() << ". Exiting!" << endl;
        return EXIT_FAILURE;
    }

    // read input fasta file
    if(read_fasta(in_data.fasta_file) != 0) {
        cout << "Error reading " << in_data.fasta_file.path << " file. Exiting!"
             << endl;
        return EXIT_FAILURE;
    } else if(in_data.fasta_file.seq_names.size() < 3) {
        cout << "At least three sequences required. Exiting!" << endl;
        return EXIT_FAILURE;
    }

    if(in_data.out_file.empty()) {  // if no output is specified save in current
                                    // dir in PHYLIP format
        in_data.out_file =
            boost::filesystem::path(in_data.fasta_file.path).stem().string() +
            ".phy";
    } else if(boost::filesystem::extension(in_data.out_file) != ".phy" &&
              boost::filesystem::extension(in_data.out_file) != ".fasta") {
        cout << "Format for output file is not valid. Exiting!" << endl;
        return EXIT_FAILURE;
    }

    return ref_indel_alignment(in_data);
}
