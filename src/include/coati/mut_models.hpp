/*
# Copyright (c) 2020 Juan J. Garcia Mesa <jgarc111@asu.edu>
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

#ifndef MUT_MODELS_H
#define MUT_MODELS_H

#include <fst/fstlib.h>
#include <Eigen/Dense>
#include <iostream>
#include <unordered_map>
#include <map>
#include <vector>
#include <algorithm>
#include <unsupported/Eigen/MatrixFunctions>
#include <coati/utils.hpp>

typedef Eigen::Matrix<double, 64, 64>Matrix64f;
typedef Eigen::Matrix<double, 64, 1>Vector64f;

using namespace fst;
using namespace std;

void nuc2pos(VectorFst<StdArc>& n2p);
void marg_mut(VectorFst<StdArc>& mut_fst, VectorFst<StdArc> marg_pos);
void toycoati(VectorFst<StdArc>& mut_fst);
void toy_marg(VectorFst<StdArc>& mut_fst);
void dna_mut(VectorFst<StdArc>& mut_fst);
void ecm_p(Matrix64f& P);
void ecm(VectorFst<StdArc>& mut_fst);
void ecm_marginal(VectorFst<StdArc>& mut_fst);
bool syn(cod c1, cod c2);
void nts_ntv(cod c1, cod c2, int& nts, int& ntv);
double k(cod c1, cod c2, int model=0);

#endif
