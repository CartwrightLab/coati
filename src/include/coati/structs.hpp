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

#ifndef STRUCTS_HPP
#define STRUCTS_HPP

#include <fst/fstlib.h>

#include <filesystem>

#include "matrix.hpp"

namespace coati {

using VectorFstStdArc = fst::VectorFst<fst::StdArc>;

struct gap_t {
    std::size_t len{1};                 /*!< unit size of gaps */
    float_t open{0.001};                /*!< gap opening score */
    float_t extend{1.0f - 1.0f / 6.0f}; /*!< gap extension score */

    gap_t() = default;
    explicit gap_t(std::size_t l, float_t o = 0.001,
                   float_t e = 1.0f - 1.0f / 6.0f)
        : len{l}, open{o}, extend{e} {}
};

// extracts extension and filename from both file.foo and ext:file.foo
struct file_type_t {
    std::string path;
    std::string type_ext;
};

struct data_t {
    std::filesystem::path path;        /*!< path to input file */
    std::vector<std::string> names;    /*!< names of fasta sequences */
    std::vector<std::string> seqs;     /*!< fasta sequences */
    float_t weight{0.f};               /*!< alignment weight */
    std::vector<VectorFstStdArc> fsts; /*!< sequences as FSTs */
    file_type_t out_file;              /*!< path to alignment output file */

    data_t() = default;
    explicit data_t(std::filesystem::path p, std::vector<std::string> n = {},
                    std::vector<std::string> s = {}, float_t w = 0.f,
                    std::vector<VectorFstStdArc> f = {},
                    std::filesystem::path o = {})
        : path{std::move(p)},
          names{std::move(n)},
          seqs{std::move(s)},
          weight{std::move(w)},
          fsts{std::move(f)},
          out_file{std::move(o)} {}

    /** \brief Return number of names/sequences */
    std::size_t size() {
        if(names.size() != seqs.size()) {
            throw std::invalid_argument(
                "Different number of sequences and names.");
        }
        return names.size();
    }

    /** \brief Return number of names/sequences */
    [[nodiscard]] size_t size() const {
        if(names.size() != seqs.size()) {
            throw std::invalid_argument(
                "Different number of sequences and names.");
        }
        return names.size();
    }

    /** \brief Return length of sequence on position index */
    size_t len(size_t index) { return seqs[index].length(); }
};

struct alignment_t {
    coati::data_t data;           /*!< sequences */
    std::string model{"m-coati"}; /*!< substitution model */
    float_t br_len{0.0133};       /*!< branch length */
    float_t omega{0.2};           /*!< nonsynonymous-synonymous bias */
    std::vector<float_t> pi{0.308, 0.185, 0.199,
                            0.308}; /*!< nucleotide frequencies */
    std::string tree{""};           /*!< path to input newick tree file */
    std::string ref{""};            /*!< name of reference sequence */
    std::string rate{""}; /*!< path to csv input substitution matrix file */
    gap_t gap;            /*!< gap struct */
    std::vector<float_t> sigma{0.f, 0.f, 0.f,
                               0.f, 0.f, 0.f}; /*!< GTR sigma parameters */
    Matrixf subst_matrix;                      /*!< substitution matrix */
    VectorFstStdArc subst_fst;                 /*!< substitution FST */
    std::filesystem::path output;      /*!< path to alignment output file */
    std::filesystem::path weight_file; /*!< file to output alignment weight */
    bool score{false}; /*!< if true an input alignment is scored */

    /** \brief Return true if model selected is marginal (m-coati or m-ecm) */
    bool is_marginal() {
        return (model.compare("m-coati") == 0 || model.compare("m-ecm") == 0);
    }
};

struct args_t {
    coati::alignment_t aln;    /*!< input data and alignment parameters */
    float_t temperature{1.0f}; /*!< temperature parameter for sampling */
    size_t sample_size{1};     /*!< sampling sample size */
};

}  // namespace coati
#endif
