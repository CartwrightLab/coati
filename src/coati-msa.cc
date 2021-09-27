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

#include <CLI11.hpp>
#include <coati/align.hpp>
#include <coati/fasta.hpp>
#include <regex>

int main(int argc, char* argv[]) {
    coati::utils::args_t args;

    // Parse command line options
    CLI::App msa;
    coati::utils::set_cli_options(msa, args, "msa");
    CLI11_PARSE(msa, argc, argv);

    // if no output is specified save in current dir in PHYLIP format
    if(args.output.empty()) {
        args.output =
            std::filesystem::path(args.fasta.path).stem().string() + ".phy";
    } else {  // check format is valid (phylip/fasta)
        const std::string extension =
            std::filesystem::path(args.output).extension();
        const std::regex valid_ext("[.phy,.fasta.fa]");
        if(!std::regex_match(extension, valid_ext)) {
            throw std::invalid_argument(
                "Output file format is invalid. Phylip and fasta files "
                "supported.");
        }
    }

    // read input fasta file
    args.fasta = coati::read_fasta(args.fasta.path.string());

    if(args.fasta.size() < 3) {
        throw std::invalid_argument(
            "At least three sequences required. Exiting!");
    }

    return ref_indel_alignment(args);
}
