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

#ifndef MUTATION_FST_HPP
#define MUTATION_FST_HPP

#include <string_view>

#include "matrix.hpp"
#include "structs.hpp"
#include "utils.hpp"

namespace coati {

// Create Muse and Gaut codon model FST.
VectorFstStdArc mg94(float br_len, float omega,
                     const std::vector<coati::float_t>& pi,
                     const std::vector<coati::float_t>& sigma = {0, 0, 0, 0, 0,
                                                                 0});
// Create dna marginal Muse and Gaut codon model FST.
VectorFstStdArc dna(float br_len, float omega,
                    const std::vector<coati::float_t>& pi);
// Create affine gap indel model FST.
VectorFstStdArc indel(float gap_open, float gap_extend,
                      const std::vector<float>& pi, float bc_error);
// Add arc to FST.
void add_arc(VectorFstStdArc& fst, int src, int dest, int ilabel = 0,
             int olabel = 0, float score = 1.0);
// Create FSA (acceptor) from a sequence.
bool acceptor(const std::string_view content, VectorFstStdArc& accept);
// Optimize FST: remove epsilons, determinize, and minimize.
VectorFstStdArc optimize(VectorFstStdArc& fst_raw);
}  // namespace coati
#endif
