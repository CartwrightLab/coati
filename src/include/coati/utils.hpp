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

#ifndef UTILS_HPP
#define UTILS_HPP

#include <fst/fstlib.h>

#include <CLI11.hpp>
#include <Eigen/Dense>
#include <boost/algorithm/string.hpp>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <regex>
#include <string>
#include <unsupported/Eigen/MatrixFunctions>
#include <vector>

#include "dna_syms.hpp"
#include "fasta.hpp"
#include "json.hpp"
#include "matrix.hpp"
#include "mg94q.tcc"
#include "mutation_coati.hpp"
#include "mutation_ecm.hpp"
#include "mutation_fst.hpp"
#include "phylip.hpp"
#include "structs.hpp"

/* Table for converting a nucleotide character to 2-bit encoding */
const uint8_t nt4_table[256] = {
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 0, 4, 1, 4, 4, 4, 2,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4};

/* Table for looking up a codon ECM group */
const uint8_t amino_group_table[64] = {
    75, 78, 75, 78, 84, 84, 84, 84, 82, 83, 82, 83, 73, 73, 77, 73,
    81, 72, 81, 72, 80, 80, 80, 80, 82, 82, 82, 82, 76, 76, 76, 76,
    69, 68, 69, 68, 65, 65, 65, 65, 71, 71, 71, 71, 86, 86, 86, 86,
    42, 89, 42, 89, 83, 83, 83, 83, 42, 67, 87, 67, 76, 70, 76, 70};

namespace coati::utils {
using VectorFstStdArc = fst::VectorFst<fst::StdArc>;

using sequence_pair_t = std::vector<std::basic_string<unsigned char>>;

coati::Matrixf parse_matrix_csv(const std::string& file);

enum struct Command { ALIGNPAIR, MSA, SAMPLE, FORMAT };

int cod_distance(uint8_t cod1, uint8_t cod2);
int cod_int(const std::string& codon);
void set_cli_options(CLI::App& app, coati::args_t& args,
                     const coati::utils::Command& command);
void set_cli_options_format(CLI::App& app, coati::args_t& args,
                            const coati::utils::Command& command);
sequence_pair_t marginal_seq_encoding(const std::string& anc,
                                      const std::string& des);
void set_subst(alignment_t& aln);

// returns {.ext, file.foo}
// trims whitespace as well
file_type_t extract_file_type(std::string path);

data_t read_input(alignment_t& aln);
bool write_output(data_t& data, const VectorFstStdArc& aln = {});
void fst_to_seqs(coati::data_t& data, const VectorFstStdArc& aln);

// calculate log(1+exp(x))
// https://cran.r-project.org/web/packages/Rmpfr/vignettes/log1mexp-note.pdf
static constexpr inline float_t log1p_exp(float_t x) {
    if(x <= -37.0f) {
        return std::exp(x);
    }
    if(x <= 18.0f) {
        return std::log1p(std::exp(x));
    }
    if(x <= 33.3f) {
        return x + std::exp(-x);
    }
    return x;
}
// calculate log(exp(a)+exp(b))
// Let x = max(a,b)
// Let y = -abs(a-b)
//  log(exp(a)+exp(b)) = x+log(1+exp(y))
static constexpr inline float_t log_sum_exp(float_t a, float_t b) {
    float_t x = std::max(a, b);
    float_t y = -std::fabs(a - b);
    return x + log1p_exp(y);
}

}  // namespace coati::utils
#endif
