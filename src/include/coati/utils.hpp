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

using VectorFstStdArc = fst::VectorFst<fst::StdArc>;

using sequence_pair_t = std::vector<std::basic_string<unsigned char>>;

namespace coati::utils {
struct gap_t {
    std::size_t len{1};
    float_t open{0.001};
    float_t extend{1.f - 1.f / 6.f};

    gap_t() = default;
    gap_t(std::size_t l, float_t o = 0.001, float_t e = 1.f - 1.f / 6.f)
        : len{l}, open{o}, extend{e} {}
};

struct args_t {
    coati::fasta_t fasta;
    std::string model{"m-coati"}, weight_file{""};
    std::filesystem::path output;
    bool score{false};
    std::string tree{""}, ref{""}, rate{""};
    gap_t gap;
    float_t br_len{0.0133};
    float_t omega{0.2};
    std::vector<float_t> pi{0.308, 0.185, 0.199, 0.308};

    args_t() = default;
    args_t(const std::string& f, const std::vector<std::string>& n,
           const std::vector<std::string>& s, std::string m = "m-coati",
           std::string weight = "", std::filesystem::path out = "",
           bool sc = false, std::string tr = "", std::string re = "",
           std::string ra = "", size_t gl = 1, float_t go = 0.001,
           float_t ge = 1.f - 1.f / 6.f, float_t br = 0.0133, float_t w = 0.2f,
           std::vector<float> p = {0.308, 0.185, 0.199, 0.308})
        : fasta{coati::fasta_t(f, n, s)},
          model{std::move(m)},
          weight_file{std::move(weight)},
          output{std::move(out)},
          score{sc},
          tree{std::move(tr)},
          ref{std::move(re)},
          rate{std::move(ra)},
          gap{gap_t(gl, go, ge)},
          br_len{br},
          omega{w},
          pi{std::move(p)} {}
};
}  // namespace coati::utils

struct alignment_t {
    coati::fasta_t fasta;
    float_t weight{0.0};
    std::filesystem::path weight_file;
    std::string model{""};

    alignment_t() = default;
    // NOLINTNEXTLINE(misc-unused-parameters)
    alignment_t(const std::filesystem::path& f,
                const std::vector<std::string>& n,
                // NOLINTNEXTLINE(misc-unused-parameters)
                const std::vector<std::string>& s, float_t w,
                std::filesystem::path w_f, std::string m)
        : fasta{coati::fasta_t(f, n, s)},
          weight{w},
          weight_file{std::move(w_f)},
          model{std::move(m)} {}
};

bool write_phylip(coati::fasta_t& fasta);
bool write_phylip(VectorFstStdArc& aln, coati::fasta_t& fasta);
int cod_distance(uint8_t cod1, uint8_t cod2);
int cod_int(const std::string& codon);
coati::Matrix<coati::float_t> parse_matrix_csv(const std::string& file);
void set_cli_options(CLI::App& app, coati::utils::args_t& in_data,
                     const std::string& command);
sequence_pair_t marginal_seq_encoding(const std::string& anc,
                                      const std::string& des);
#endif
