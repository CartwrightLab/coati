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

#ifndef MUTATION_COATI_HPP
#define MUTATION_COATI_HPP

#include <algorithm>
#include <cmath>
#include <limits>
#include <unsupported/Eigen/MatrixFunctions>
#include <vector>

#include "matrix.hpp"
#include "mg94p.tcc"
#include "utils.hpp"

namespace coati {
// Create Muse \& Gaut (1994) substitution matrix.
coati::Matrixf mg94_p(coati::float_t br_len, coati::float_t omega,
                      const std::vector<coati::float_t>& nuc_freqs,
                      const std::vector<coati::float_t>& sigma = {0, 0, 0, 0, 0,
                                                                  0});
// Create marginal 192x4 substitution P matrix.
coati::Matrixf marginal_p(const coati::Matrixf& P,
                          const std::vector<coati::float_t>& pi,
                          const coati::AmbiguousNucs amb);
// Probabilities for ambiguous nucs in marginal model using averages.
void ambiguous_avg_p(coati::Matrixf& p);
// Probabilities for ambiguous nucs in marginal model taking best prob.
void ambiguous_best_p(coati::Matrixf& p);
// Create GTR substitution model matrix.
coati::Matrixf gtr_q(const std::vector<coati::float_t>& nuc_freqs,
                     const std::vector<coati::float_t>& sigma);
}  // namespace coati
#endif
