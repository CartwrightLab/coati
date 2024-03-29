/*
# Copyright (c) 2020-2022 Juan J. Garcia Mesa <juanjosegarciamesa@gmail.com>
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

#include <fst/fst-decl.h>

#include <CLI11.hpp>
#include <string>
#include <string_view>
#include <vector>

#include "structs.hpp"

/**
 * @brief Table for converting a nucleotide character to 4-bit encoding.
 *
 * Following IUPAC nucleotide code
 * | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 | 15 |
 * | A | C | G |T/U| R | Y | M | K | S | W | B  | D  | H  | V  | N  | -  |
 *
 * R: purine       A or G
 * Y: pyrimidine   C or T/U
 * M: amino group  A or C
 * K: keto group   G or T/U
 * S: strong inter C or G
 * W: weak interac A or T/U
 * B: not A        C or G or T/U
 * D: not C        A or G or T/U
 * H: not G        A or C or T/U
 * V: not T        A or C or G
 * N: any          A or C or G or T/U
 */
const uint8_t nt16_table[128] = {
    16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16,
    16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16,
    16, 16, 16, 16, 16, 16, 16, 15, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16,
    16, 16, 16, 16, 16, 16, 16, 16, 0,  10, 1,  11, 16, 16, 2,  12, 16, 16, 7,
    16, 6,  14, 16, 16, 16, 4,  8,  3,  3,  13, 9,  16, 5,  16, 16, 16, 16, 16,
    16, 16, 0,  10, 1,  11, 16, 16, 2,  12, 16, 16, 7,  16, 6,  14, 16, 16, 16,
    4,  8,  3,  3,  13, 9,  16, 5,  16, 16, 16, 16, 16};

/**
 * @brief Table for looking up a codon ECM group.
 */
const uint8_t amino_group[61] = {
    75, 78, 75, 78, 84, 84, 84, 84, 82, 83, 82, 83, 73, 73, 77, 73,
    81, 72, 81, 72, 80, 80, 80, 80, 82, 82, 82, 82, 76, 76, 76, 76,
    69, 68, 69, 68, 65, 65, 65, 65, 71, 71, 71, 71, 86, 86, 86, 86,
    89, 89, 83, 83, 83, 83, 67, 87, 67, 76, 70, 76, 70};

namespace coati::utils {

using sequence_pair_t = std::vector<std::basic_string<unsigned char>>;

// Hamming distance between two codons.
int cod_distance(uint8_t cod1, uint8_t cod2);
// Get a codon's position in the codon list (AAA->0, AAC->1, ..., TTT->63).
int cod_int(const std::string_view codon);
// Setup command line options for coati-alignpair.
void set_options_alignpair(CLI::App& app, coati::args_t& args);
// Setup command line options for coati-msa.
void set_options_msa(CLI::App& app, coati::args_t& args);
// Setup command line options for coati-sample.
void set_options_sample(CLI::App& app, coati::args_t& args);
// Setup command line options for coati format.
void set_options_format(CLI::App& app, coati::args_t& args);
// Encode two sequences as vector<unsigned char>.
sequence_pair_t marginal_seq_encoding(const std::string_view anc,
                                      const std::string_view des);
// Set substitution matrix or FST according to model.
void set_subst(alignment_t& aln);
// Extract file type from path. Returns {.ext, file.foo}.
// trims whitespace as well
file_type_t extract_file_type(std::string path);
// Convert alignment FST to std::string sequences.
void fst_to_seqs(coati::data_t& data, const VectorFstStdArc& aln);
// Get nucleotide from codon list without stop codons.
uint8_t get_nuc(uint8_t cod, int pos);
// Reorder pair of input sequences so that reference is at position zero
void order_ref(coati::alignment_t& aln);
// Validate input sequences
void process_marginal(coati::alignment_t& aln);
// Validate input pairwise alignment for scoring
std::string process_alignment(coati::alignment_t& aln);
// Trim end stop codons
void trim_end_stops(coati::data_t& data);
// Restore end stop codons
void restore_end_stops(coati::data_t& data, const coati::gap_t& gap);
// Read and validate input sequences
void process_triplet(coati::alignment_t& aln);
// Convert codon index from 64 codon table to 61 (no stop codons)
int cod64_to_61(int cod);
// Convert codon index from 61 (no stop codons) codon table to 64
int cod61_to_64(int cod);

// calculate log(1+exp(x))
// https://cran.r-project.org/web/packages/Rmpfr/vignettes/log1mexp-note.pdf
// https://cs.stackexchange.com/questions/110798/numerically-stable-log1pexp-calculation
inline double log1p_exp(double x) {
    using namespace std;
    if(x <= -37.0) {  // log(epsilon)
        return exp(x);
    }
    if(x <= 18.0) {  // -log(2 * .Machine$double.eps) / 2
        return log1p(exp(x));
    }
    if(x <= 33.3) {
        return x + exp(-x);
    }
    return x;
}

inline float log1p_exp(float x) {
    using namespace std;
    if(x <= -16.0f) {  // log(epsilon)
        return expf(x);
    }
    if(x <= 8.0f) {  // -log(2 * epsilon) / 2
        return log1pf(expf(x));
    }
    if(x <= 14.5f) {
        return x + expf(-x);
    }
    return x;
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

inline float_t log_sum_exp(float_t a, float_t b, float_t c) {
    return log_sum_exp(log_sum_exp(a, b), c);
}

}  // namespace coati::utils
#endif
