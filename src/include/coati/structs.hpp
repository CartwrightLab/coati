/*
# Copyright (c) 2021-2022 Juan J. Garcia Mesa <juanjosegarciamesa@gmail.com>
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

#include <filesystem>

#include "data.hpp"
#include "matrix.hpp"

namespace coati {

/**
 * @brief Stores gap unit length and open/extend score cost.
 *
 */
struct gap_t {
   public:
    std::size_t len{1};                 /*!< unit length of gaps */
    float_t open{0.001};                /*!< gap opening score */
    float_t extend{1.0f - 1.0f / 6.0f}; /*!< gap extension score */

    gap_t() = default;
    explicit gap_t(std::size_t l, float_t o = 0.001,
                   float_t e = 1.0f - 1.0f / 6.0f)
        : len{l}, open{o}, extend{e} {}
};

/**
 * @brief Extracts extension and filename from both file.foo and ext:file.foo
 *
 */
struct file_type_t {
   public:
    std::string path;
    std::string type_ext;
};

enum struct AmbiguousNucs { SUM, BEST };
enum struct MarginalSubst { SUM, MAX };

/**
 * @brief Stores input and model parameters from an alignment.
 *
 */
class alignment_t {
   public:
    coati::data_t data;          /*!< sequences */
    std::string model{"mar-mg"}; /*!< substitution model */
    float_t br_len{0.0133};      /*!< branch length */
    float_t omega{0.2};          /*!< nonsynonymous-synonymous bias */
    std::vector<float_t> pi{0.308, 0.185, 0.199,
                            0.308}; /*!< nucleotide frequencies */
    std::string tree{""};           /*!< path to input newick tree file */
    std::string refs{""};           /*!< name of reference sequence */
    bool rev{false};                /*!< use 2nd seq as reference */
    std::string rate{""};           /*!< path to csv input subst matrix file */
    gap_t gap;                      /*!< gap struct */
    std::vector<float_t> sigma{0.f, 0.f, 0.f,
                               0.f, 0.f, 0.f}; /*!< GTR sigma parameters */
    Matrixf subst_matrix;                      /*!< substitution matrix */
    VectorFstStdArc subst_fst;                 /*!< substitution FST */
    std::filesystem::path output; /*!< path to alignment output file */
    bool score{false};            /*!< if true an input alignment is scored */
    AmbiguousNucs amb = AmbiguousNucs::SUM;
    MarginalSubst sub = MarginalSubst::SUM;
    float_t bc_error{0.0001}; /*!< base calling error rate */

    /** \brief Return true if model selected is marginal (marginal or m-ecm) */
    bool is_marginal() {
        return (model == "mar-mg" || model == "mar-ecm" || !rate.empty());
    }

    /** \brief Return sequence at position index. */
    std::string& seq(size_t index) { return data.seqs[index]; }

    /** \brief Return name of sequence at position index. */
    std::string& name(size_t index) { return data.names[index]; }
};

/**
 * @brief Stores parameters for running coati format.
 *
 */
struct format_t {
   public:
    bool preserve_phase{false}; /*!< preserve phase for downstream analyses */
    std::string padding{"?"};   /*!< padding char to format preserve phase */
    std::vector<std::string> names{}; /*!< names of seqs to extract */
    std::vector<size_t> pos{};        /*!< position of seqs to extract */
};

/**
 * @brief Stores parameters for running coati sample.
 *
 */
struct sample_t {
    float_t temperature{1.0f}; /*!< temperature parameter for sampling */
    size_t sample_size{1};     /*!< sampling sample size */
    std::vector<std::string> seeds{{""}}; /*!< seeds for sampling */
};

/**
 * @brief Stores all input values for any coati command.
 *
 */
class args_t {
   public:
    coati::alignment_t aln; /*!< input data and alignment parameters */
    coati::sample_t sample; /*!< coati sample arguments */
    coati::format_t format; /*!< coati format arguments */
};

}  // namespace coati
#endif
