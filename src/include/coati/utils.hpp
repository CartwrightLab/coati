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
#include "matrix.hpp"
#include "mg94q.tcc"
#include "mutation_coati.hpp"
#include "mutation_ecm.hpp"
#include "mutation_fst.hpp"

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

struct gap_t {
    std::size_t len{1};                 /*!< unit size of gaps */
    float_t open{0.001};                /*!< gap opening score */
    float_t extend{1.0f - 1.0f / 6.0f}; /*!< gap extension score */

    gap_t() = default;
    explicit gap_t(std::size_t l, float_t o = 0.001,
                   float_t e = 1.0f - 1.0f / 6.0f)
        : len{l}, open{o}, extend{e} {}
};

struct args_t {
    coati::fasta_t fasta;         /*!< fasta struct */
    std::string model{"m-coati"}; /*!< substitution model */
    std::string weight_file{""};  /*!< file to output alignment weight */
    std::filesystem::path output; /*!< path to alignment output file */
    bool score{false};            /*!< if true an input alignment is scored */
    std::string tree{""};         /*!< path to input newick tree file */
    std::string ref{""};          /*!< name of reference sequence */
    std::string rate{""};   /*!< path to csv input substitution matrix file */
    gap_t gap;              /*!< gap struct */
    float_t br_len{0.0133}; /*!< branch length */
    float_t omega{0.2};     /*!< nonsynonymous-synonymous bias */
    std::vector<float_t> pi{0.308, 0.185, 0.199,
                            0.308}; /*!< nucleotide frequencies */
    float_t temperature{1.0f};
    size_t sample_size{1};
};

coati::Matrixf parse_matrix_csv(const std::string& file);

struct alignment_t {
    coati::fasta_t fasta;              /*!< fasta struct */
    float_t weight{0.0};               /*!< alignment weight */
    std::filesystem::path weight_file; /*!< file to output alignment weight */
    std::string model{""};             /*!< substitution model */
    Matrixf subst_matrix;              /*!< substitution matrix */
    VectorFstStdArc subst_fst;         /*!< substitution FST */
    std::vector<VectorFstStdArc> seqs = {}; /*!< sequences as FSTs */

    alignment_t() = default;
    // NOLINTNEXTLINE(misc-unused-parameters)
    alignment_t(const std::filesystem::path& f,
                const std::vector<std::string>& n,
                // NOLINTNEXTLINE(misc-unused-parameters)
                const std::vector<std::string>& s, float_t w,
                std::filesystem::path w_f, std::string m, Matrixf p,
                VectorFstStdArc fst, std::vector<VectorFstStdArc> ss)
        : fasta{coati::fasta_t(f, n, s)},
          weight{w},
          weight_file{std::move(w_f)},
          model{std::move(m)},
          subst_matrix{std::move(p)},
          subst_fst{std::move(fst)},
          seqs{std::move(ss)} {}

    /** \brief Return true if model selected is marginal (m-coati or m-ecm) */
    bool is_marginal() {
        return (model.compare("m-coati") == 0 || model.compare("m-ecm") == 0);
    }
};

enum struct Command { ALIGNPAIR, MSA, SAMPLE, FORMAT };

int cod_distance(uint8_t cod1, uint8_t cod2);
int cod_int(const std::string& codon);
void set_cli_options(CLI::App& app, coati::utils::args_t& args,
                     const coati::utils::Command& command);
sequence_pair_t marginal_seq_encoding(const std::string& anc,
                                      const std::string& des);
void set_subst(args_t& args, alignment_t& aln);

// extracts extension and filename from both file.foo and ext:file.foo
struct file_type_t {
    std::string path;
    std::string type_ext;
};

// returns {.ext, file.foo}
// trims whitespace as well
file_type_t extract_file_type(std::string path);

// calculate log(1+exp(x))
// https://cran.r-project.org/web/packages/Rmpfr/vignettes/log1mexp-note.pdf
inline float_t log1p_exp(float_t x) {
    if(x <= -37.0f) {
        return std::exp(x);
    } else if(x <= 18.0f) {
        return std::log1p(std::exp(x));
    } else if(x <= 33.3f) {
        return x + std::exp(-x);
    } else {
        return x;
    }
}
// calculate log(exp(a)+exp(b))
// Let x = max(a,b)
// Let y = -abs(a-b)
//  log(exp(a)+exp(b)) = x+log(1+exp(y))
inline float_t log_sum_exp(float_t a, float_t b) {
    float_t x = std::max(a, b);
    float_t y = -std::fabs(a - b);
    return x + log1p_exp(y);
}

}  // namespace coati::utils
#endif
