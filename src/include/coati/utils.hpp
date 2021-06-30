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

#include <Eigen/Dense>
#include <boost/algorithm/string.hpp>
#include <coati/dna_syms.hpp>

#include <fst/fstlib.h>
#include <fstream>
#include <iostream>
#include <filesystem>
#include <string>
#include <vector>

/* Table for converting a nucleotide character to 2-bit encoding and for
        looking up coding amino acid based on codon index (AAA:0, AAC:1, ...,
   TTT:63) */
const uint8_t nt4_table[256] = {
    75, 78, 75, 78, 84, 84, 84, 84, 82, 83, 82, 83, 73, 73, 77, 73, 81, 72, 81,
    72, 80, 80, 80, 80, 82, 82, 82, 82, 76, 76, 76, 76, 69, 68, 69, 68, 65, 65,
    65, 65, 71, 71, 71, 71, 86, 86, 86, 86, 42, 89, 42, 89, 83, 83, 83, 83, 42,
    67, 87, 67, 76, 70, 76, 70, 4,  0,  4,  1,  4,  4,  4,  2,  4,  4,  4,  4,
    4,  4,  4,  4,  4,  4,  4,  4,  3,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,
    4,  4,  0,  4,  1,  4,  4,  4,  2,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,
    4,  4,  3,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,
    4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,
    4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,
    4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,
    4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,
    4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,
    4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,
    4,  4,  4,  4,  4,  4,  4,  4,  4};

typedef Eigen::Matrix<double, 64, 64> Matrix64f;
typedef Eigen::Matrix<double, 5, 1> Vector5d;

struct fasta_t {
    std::filesystem::path path;
    std::vector<std::string> seq_names, seq_data;
    // fasta_t() : path{""}, seq_names{}, seq_data{} {}
    fasta_t(std::string f = "", std::vector<std::string> n = {},
            const std::vector<std::string> d = {})
        : path{std::move(f)}, seq_names{std::move(n)}, seq_data{std::move(d)} {}
};

struct input_t {
    std::string mut_model, weight_file, out_file, rate, tree, ref;
    bool score;
    double br_len;
    fasta_t fasta_file;
};

struct alignment_t {
    std::string weight_file, model;
    float weight;
    fasta_t f;
    alignment_t() : f{fasta_t()}, weight{0.0}, weight_file{""}, model{""} {}
    alignment_t(const std::string& f_file, const std::vector<std::string>& n,
                const std::vector<std::string>& d, float w, const std::string& w_f,
                const std::string& m)
        : f{fasta_t(f_file, n, d)}, weight{w}, weight_file{w_f}, model{m} {}
};

int read_fasta(fasta_t& fasta_file, std::vector<fst::VectorFst<fst::StdArc>>& fsts);
int read_fasta(fasta_t& fasta_file);
void add_arc(fst::VectorFst<fst::StdArc>& fst, int src, int dest, int ilabel = 0,
             int olabel = 0, float weight = 1.0);
fst::VectorFst<fst::StdArc> optimize(fst::VectorFst<fst::StdArc> fst);
int write_fasta(fasta_t& fasta_file);
int write_fasta(fst::VectorFst<fst::StdArc>& aln, fasta_t& fasta_file);
int write_phylip(fasta_t& fasta_file);
int write_phylip(fst::VectorFst<fst::StdArc>& aln, fasta_t& fasta_file);
bool acceptor(std::string content, fst::VectorFst<fst::StdArc>& accept);
int cod_distance(uint8_t cod1, uint8_t cod2);
int cod_int(std::string codon);
int parse_matrix_csv(std::string file, Matrix64f& P, double& br_len);

#endif
