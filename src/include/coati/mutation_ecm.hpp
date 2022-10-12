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

#ifndef MUTATION_ECM_HPP
#define MUTATION_ECM_HPP

#include "ecm.tcc"
#include "matrix.hpp"
#include "mutation_fst.hpp"

namespace coati {
using VectorFstStdArc = fst::VectorFst<fst::StdArc>;

// Calculate number of transitions and transversion between two codons.
void nts_ntv(uint8_t c1, uint8_t c2, int& nts, int& ntv);
// Transition-transverison bias function.
float k(uint8_t c1, uint8_t c2, int model = 0, float kappa = 2.5);
// Create Empirical Codon Model substitution P matrix.
coati::Matrixf ecm_p(float br_len, float omega);
// Create Empirical Codon Model (Kosiol et al. 2007) FST.
VectorFstStdArc ecm(float br_len, float omega);
}  // namespace coati
#endif
