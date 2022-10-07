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

#ifndef ALIGN_PAIR_HPP
#define ALIGN_PAIR_HPP

#include <cmath>
#include <limits>
#include <random.hpp>
#include <string_view>

#include "semiring.hpp"
#include "structs.hpp"
#include "utils.hpp"

namespace coati {

using random_t = fragmites::random::Random;
using base_t = unsigned char;
using seq_view_t = std::basic_string_view<base_t>;
using stationary_vector_t = std::vector<float_t>;

/**
 * @brief Stores matrices for all calculations in dynamic programming alignment.
 *
 */

class align_pair_work_t {
   public:
    Matrixf mch;     /*!< match state */
    Matrixf del;     /*!< deletion state */
    Matrixf ins;     /*!< insertion state */
    Matrixf mch_mch; /*!< match to match state */
    Matrixf mch_del; /*!< match to del state */
    Matrixf mch_ins; /*!< match to ins state */
    Matrixf del_mch; /*!< del to mch state */
    Matrixf del_del; /*!< del to del state */
    Matrixf ins_mch; /*!< ins to mch state */
    Matrixf ins_del; /*!< ins to del state */
    Matrixf ins_ins; /*!< ins to ins state */

    float_t subst{0.0f};   /*!< substitution tmp value*/
    float_t mch2mch{0.0f}; /*!< match to match tmp value */
    float_t mch2del{0.0f}; /*!< match to del tmp value */
    float_t mch2ins{0.0f}; /*!< match to ins tmp value */
    float_t del2mch{0.0f}; /*!< del to mch tmp value */
    float_t del2del{0.0f}; /*!< del to del tmp value */
    float_t ins2mch{0.0f}; /*!< ins to mch tmp value */
    float_t ins2del{0.0f}; /*!< ins to del tmp value */
    float_t ins2ins{0.0f}; /*!< ins to ins tmp value */

    /** \brief Resize matrices and initialize to given value
     *
     * @param[in] len_a std::size_t number of rows.
     * @param[in] len_b std::size_t number of columns.
     * @param[in] val coati::float_t value.
     */
    void resize(size_t len_a, size_t len_b, float_t val) {
        mch.resize(len_a, len_b, val);      // match matrix
        del.resize(len_a, len_b, val);      // deletion matrix
        ins.resize(len_a, len_b, val);      // insertion matrix
        mch_mch.resize(len_a, len_b, val);  // match to match matrix
        mch_del.resize(len_a, len_b, val);  // match to del matrix
        mch_ins.resize(len_a, len_b, val);  // match to ins matrix
        del_mch.resize(len_a, len_b, val);  // del to match matrix
        del_del.resize(len_a, len_b, val);  // del to del matrix
        ins_mch.resize(len_a, len_b, val);  // ins to match matrix
        ins_ins.resize(len_a, len_b, val);  // ins to ins matrix
        ins_del.resize(len_a, len_b, val);  // ins to del matrix
    }

    /** \brief Save edge values between mch, del, and ins states
     *
     * @param[in] i std::size_t row index.
     * @param[in] j std::size_t column index.
     */
    void save_values(size_t i, size_t j) {
        mch_mch(i, j) = mch2mch;
        mch_del(i, j) = mch2del;
        mch_ins(i, j) = mch2ins;
        del_mch(i, j) = del2mch;
        del_del(i, j) = del2del;
        ins_mch(i, j) = ins2mch;
        ins_del(i, j) = ins2del;
        ins_ins(i, j) = ins2ins;
    }

    /** \brief Initialize margins in edge matrices del to del and ins to ins
     *
     */
    void init_margins() {
        del_del = del;
        ins_ins = ins;
    }
};

/**
 * @brief Stores matrices and temporary values in dynamic programming alignment.
 *
 */
class align_pair_work_mem_t {
   public:
    Matrixf mch;           /*!< match state */
    Matrixf del;           /*!< deletion state */
    Matrixf ins;           /*!< insertion state */
    float_t subst{0.0f};   /*!< substitution tmp value*/
    float_t mch2mch{0.0f}; /*!< match to match tmp value */
    float_t mch2del{0.0f}; /*!< match to del tmp value */
    float_t mch2ins{0.0f}; /*!< match to ins tmp value */
    float_t del2mch{0.0f}; /*!< del to mch tmp value */
    float_t del2del{0.0f}; /*!< del to del tmp value */
    float_t ins2mch{0.0f}; /*!< ins to mch tmp value */
    float_t ins2del{0.0f}; /*!< ins to del tmp value */
    float_t ins2ins{0.0f}; /*!< ins to ins tmp value */

    /** \brief Resize matrices and initialize to given value
     *
     * @param[in] len_a std::size_t number of rows.
     * @param[in] len_b std::size_t number of columns.
     * @param[in] val coati::float_t value.
     */
    void resize(size_t len_a, size_t len_b, float_t val) {
        mch.resize(len_a, len_b, val);  // match matrix
        del.resize(len_a, len_b, val);  // deletion matrix
        ins.resize(len_a, len_b, val);  // insertion matrix
    }

    void save_values(size_t i, size_t j) {}
    void init_margins() {}
};

/** \brief Max value between two coati::float_t values */
inline float_t maximum(float_t x, float_t y) { return std::max(x, y); }
/** \brief Max value between three coati::float_t values */
inline float_t maximum(float_t x, float_t y, float_t z) {
    return std::max(maximum(x, y), z);
}

template <class S, class W>
void forward_impl(W &work, const seq_view_t &a, const seq_view_t &b,
                  const alignment_t &aln);

void forward(align_pair_work_t &work, const seq_view_t &a, const seq_view_t &b,
             const alignment_t &aln);
void forward_mem(align_pair_work_mem_t &work, const seq_view_t &a,
                 const seq_view_t &b, const alignment_t &aln);
void viterbi(align_pair_work_t &work, const seq_view_t &a, const seq_view_t &b,
             const alignment_t &aln);
void viterbi_mem(align_pair_work_mem_t &work, const seq_view_t &a,
                 const seq_view_t &b, const alignment_t &aln);

void traceback(const align_pair_work_mem_t &work, const std::string &a,
               const std::string &b, alignment_t &aln, size_t look_back);

void sampleback(const align_pair_work_t &work, const std::string &a,
                const std::string &b, alignment_t &aln, size_t look_back,
                random_t &rand);

}  // namespace coati
#endif
