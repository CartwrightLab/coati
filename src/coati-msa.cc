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
            "fasta,f",
            po::value<std::filesystem::path>(&in_data.fasta_file.path)
                ->required(),
            "fasta file path")("model,m",
                               po::value<std::string>(&in_data.mut_model)
                                   ->default_value("m-coati"),
                               "substitution model: m-coati (default), m-ecm")(
            "weight,w", po::value<std::string>(&in_data.weight_file),
            "write alignment score to file (coming soon)")(
            "output,o", po::value<std::string>(&in_data.out_file),
            "alignment output file")(
            "tree,p", po::value<std::string>(&in_data.tree)->required(),
            "newick phylogenetic tree")(
            "ref,r", po::value<std::string>(&in_data.ref)->required(),
            "reference sequence")("no-frameshifts", "Don't allow frameshifts")(
            "gap-open,g", po::value<float>(&in_data.gapo)->default_value(0.001),
            "Gap opening score")(
            "gap-extend,e",
            po::value<float>(&in_data.gape)->default_value(1.f - (1.f / 6.f)),
            "Gap extension score")(
            "omega,w", po::value<float>(&in_data.omega)->default_value(0.2),
            "Nonsynonymous-synonymous bias");

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
            std::cout << "Usage: coati msa file.fasta tree.newick [options]"
                      << std::endl
                      << std::endl;
            std::cout << desc << std::endl;
            return EXIT_SUCCESS;
        }

        if(varm.count("no-frameshifts")) {
            in_data.frameshifts = false;
        }

        po::notify(varm);

    } catch(po::error& e) {
        std::cerr << e.what() << ". Exiting!" << std::endl;
        return EXIT_FAILURE;
    }

    // read input fasta file
    if(read_fasta(in_data.fasta_file) != 0) {
        throw std::invalid_argument("Error reading " +
                                    in_data.fasta_file.path.string() +
                                    " file. Exiting!");
    } else if(in_data.fasta_file.seq_names.size() < 3) {
        throw std::invalid_argument(
            "At least three sequences required. Exiting!");
    }

    // if no output is specified save in current dir in PHYLIP format
    if(in_data.out_file.empty()) {
        in_data.out_file =
            std::filesystem::path(in_data.fasta_file.path).stem().string() +
            ".phy";
    } else if(std::filesystem::path(in_data.out_file).extension() != ".phy" &&
              std::filesystem::path(in_data.out_file).extension() != ".fasta") {
        throw std::invalid_argument(
            "Format for output file is not valid. Exiting!");
    }

    return ref_indel_alignment(in_data);
}
