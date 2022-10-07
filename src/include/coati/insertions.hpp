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

#ifndef INSERTIONS_HPP
#define INSERTIONS_HPP

#include <Eigen/Sparse>
#include <iostream>
#include <numeric>
#include <string_view>
#include <vector>

namespace coati {
using SparseVectorInt = Eigen::SparseVector<int, Eigen::RowMajor>;

/**
 * @brief Keep track of open and close insertions for msa.
 *
 */
struct insertion_data_t {
    std::vector<std::string> sequences; /*!< sequences */
    std::vector<std::string> names;     /*!< sequence names */
    SparseVectorInt
        insertions; /*!< insertion positions and type (open/closed) */

    insertion_data_t() = default;

    insertion_data_t(const std::string& s, const std::string& n,
                     const SparseVectorInt& i)
        : sequences(1, s), names(1, n), insertions{i} {}

    insertion_data_t(std::vector<std::string> s, std::vector<std::string> n,
                     const SparseVectorInt& i)
        : sequences{std::move(s)}, names{std::move(n)}, insertions{i} {}
};

using insertion_vector = std::vector<insertion_data_t>;

bool insertion_flags(const std::string_view ref, const std::string_view seq,
                     SparseVectorInt& insertions_vector);
bool merge_indels(coati::insertion_vector& ins_data,
                  insertion_data_t& merged_data);
uint64_t add_closed_ins(coati::insertion_vector& ins_data, std::size_t pos);
bool check_all_open(coati::insertion_vector& ins_data, std::size_t pos);
std::vector<std::size_t> find_open_ins(coati::insertion_vector& ins_data,
                                       std::size_t pos);
void add_gap(coati::insertion_vector& ins_data,
             const std::vector<std::size_t>& seq_indexes, std::size_t pos);
}  // namespace coati
#endif
