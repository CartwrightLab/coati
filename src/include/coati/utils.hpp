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
const uint8_t codon_table[64] = {
    75, 78, 75, 78, 84, 84, 84, 84, 82, 83, 82, 83, 73, 73, 77, 73,
    81, 72, 81, 72, 80, 80, 80, 80, 82, 82, 82, 82, 76, 76, 76, 76,
    69, 68, 69, 68, 65, 65, 65, 65, 71, 71, 71, 71, 86, 86, 86, 86,
    42, 89, 42, 89, 83, 83, 83, 83, 42, 67, 87, 67, 76, 70, 76, 70};

using VectorFstStdArc = fst::VectorFst<fst::StdArc>;

struct fasta_t {
    std::filesystem::path path;
    std::vector<std::string> seq_names, seq_data;
    explicit fasta_t(std::string f = "", std::vector<std::string> n = {},
                     std::vector<std::string> d = {})
        : path{std::move(f)}, seq_names{std::move(n)}, seq_data{std::move(d)} {}
};

struct input_t {
    fasta_t fasta_file;
    std::string mut_model{"m-coati"}, weight_file{""}, out_file{""}, tree{""},
        ref{""}, rate{""};
    bool score{false}, frameshifts{true};
    float br_len{0.0133}, gapo{0.001}, gape{1.f - 1.f / 6.f}, omega{0.2};

    input_t() = default;
    input_t(const std::string& f, const std::vector<std::string>& n,
            const std::vector<std::string>& d, std::string model = "m-coati",
            std::string weight = "", std::string out = "", std::string tr = "",
            std::string re = "", std::string ra = "", bool sc = false,
            bool frame = true, float br = 0.0133, float g = 0.001,
            float e = 1.f - 1.f / 6.f, float w = 0.2)
        : fasta_file{fasta_t(f, n, d)},
          mut_model{std::move(model)},
          weight_file{std::move(weight)},
          out_file{std::move(out)},
          tree{std::move(tr)},
          ref{std::move(re)},
          rate{std::move(ra)},
          score{sc},
          frameshifts{frame},
          br_len{br},
          gapo{g},
          gape{e},
          omega{w} {}
};

struct alignment_t {
    fasta_t f;
    float weight{0.0};
    std::string weight_file{""}, model{""};

    alignment_t() = default;
    // NOLINTNEXTLINE(misc-unused-parameters)
    alignment_t(const std::string& f_file, const std::vector<std::string>& n,
                // NOLINTNEXTLINE(misc-unused-parameters)
                const std::vector<std::string>& d, float w, std::string w_f,
                std::string m)
        : f{fasta_t(f_file, n, d)},
          weight{w},
          weight_file{std::move(w_f)},
          model{std::move(m)} {}
};

int read_fasta(fasta_t& fasta_file, std::vector<VectorFstStdArc>& fsts);
int read_fasta(fasta_t& fasta_file);
void add_arc(VectorFstStdArc& fst, int src, int dest, int ilabel = 0,
             int olabel = 0, float weight = 1.0);
VectorFstStdArc optimize(VectorFstStdArc fst_raw);
int write_fasta(fasta_t& fasta_file);
int write_fasta(VectorFstStdArc& aln, fasta_t& fasta_file);
int write_phylip(fasta_t& fasta_file);
int write_phylip(VectorFstStdArc& aln, fasta_t& fasta_file);
bool acceptor(std::string content, VectorFstStdArc& accept);
int cod_distance(uint8_t cod1, uint8_t cod2);
int cod_int(std::string codon);
Matrix parse_matrix_csv(const std::string& file);
void read_fasta_pair(fasta_t& fasta_file, std::vector<VectorFstStdArc>& fsts,
                     bool fst);
void set_cli_options(CLI::App& app, input_t& in_data,
                     const std::string& command);
#endif
