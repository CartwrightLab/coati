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

#ifndef ALIGN_MSA_HPP
#define ALIGN_MSA_HPP

#include "align_pair.hpp"
#include "insertions.hpp"
#include "io.hpp"
#include "semiring.hpp"
#include "tree.hpp"
#include "utils.hpp"

namespace coati {
// Multiple sequence alignment using an iterative algorithm.
bool ref_indel_alignment(coati::alignment_t& input);
// Pairwise alignments of leafs with reference sequence.
void align_leafs(coati::alignment_t& input, const coati::tree::tree_t& tree,
                 std::size_t ref_pos, const std::string& ref_seq,
                 coati::insertion_vector& nodes_ins);
// Merge alignments starting from leafs until root.
void merge_alignments(std::vector<bool>& visited,
                      const coati::tree::tree_t& tree,
                      coati::insertion_vector& nodes_ins,
                      const std::vector<std::size_t>& inode_indexes);
}  // namespace coati
#endif
