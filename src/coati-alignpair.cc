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
#include <coati/align.hpp>
#include <coati/utils.hpp>

int main(int argc, char* argv[]) {
    input_t in_data;

    // Parse command line options
    CLI::App alignpair;
    set_cli_options(alignpair, in_data, "alignpair");
    CLI11_PARSE(alignpair, argc, argv);

    std::vector<VectorFstStdArc> fsts;

    if(in_data.mut_model.compare("m-coati") == 0 ||
       in_data.mut_model.compare("m-ecm") == 0) {
        read_fasta_pair(in_data.fasta_file, fsts, false);
        return mcoati(in_data);
    } else {
        read_fasta_pair(in_data.fasta_file, fsts, true);
        return fst_alignment(in_data, fsts);
    }
    return EXIT_FAILURE;
}
