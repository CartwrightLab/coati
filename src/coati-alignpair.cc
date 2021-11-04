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

#include <CLI11.hpp>
#include <coati/align_fst.hpp>
#include <coati/align_marginal.hpp>
#include <coati/fasta.hpp>
#include <coati/utils.hpp>

int main(int argc, char* argv[]) {
    coati::utils::args_t args;

    // Parse command line options
    CLI::App alignpair;
    coati::utils::set_cli_options(alignpair, args,
                                  coati::utils::Command::ALIGNPAIR);
    CLI11_PARSE(alignpair, argc, argv);

    // if no output is specified save in current dir in PHYLIP format
    if(args.output.empty()) {
        args.output = args.fasta.path.stem();
        args.output += std::filesystem::path(".phy");
    } else {  // check format is valid (phylip/fasta)
        const std::string extension = args.output.extension();
        const std::regex valid_ext("^.phy$|^.fasta$|^.fa$");
        if(!std::regex_match(extension, valid_ext)) {
            throw std::invalid_argument(
                "Output file format is invalid. Phylip and fasta files "
                "supported.");
        }
    }

    coati::utils::alignment_t aln;
    coati::utils::set_subst(args, aln);

    // subst models aligned by dynamic programming
    if(aln.is_marginal()) {
        args.fasta = coati::read_fasta(args.fasta.path.string());
        aln.fasta.path = args.output;
        aln.fasta.names = args.fasta.names;

        if(args.fasta.size() != 2) {
            throw std::invalid_argument("Exactly two sequences required.");
        }
        return coati::marg_alignment(args, aln) ? 0 : 1;
    }
    // subst models aligned by FST composition
    args.fasta = coati::read_fasta(args.fasta.path.string(), aln.seqs);
    aln.fasta.path = args.output;
    aln.fasta.names = args.fasta.names;

    if(args.fasta.names.size() != 2 || args.fasta.size() != aln.seqs.size()) {
        throw std::invalid_argument("Exactly two sequences required.");
    }
    return coati::fst_alignment(args, aln) ? 0 : 1;
}
