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

#ifndef PROFILE_ALN_H
#define PROFILE_ALN_H

#include <Eigen/Dense>
#include <coati/mut_models.hpp>

typedef Eigen::Matrix<double, 4, 1> Vector4d;
typedef Eigen::Matrix<double, 4, 3> Matrix4x3d;

Eigen::MatrixXd create_profile(string seq);
Eigen::MatrixXd create_profile(vector<string>& aln);
double transition(Matrix4x3d cod, int pos, Vector4d nuc,
                  const Eigen::Tensor<double, 3>& p);
int gotoh_profile_marginal(vector<string> seqs1, vector<string> seqs2,
                           alignment_t& aln, Matrix64f& P_m);
int backtracking_profile(Eigen::MatrixXi Bd, Eigen::MatrixXi Bp,
                         Eigen::MatrixXi Bq, vector<string> seqs1,
                         vector<string> seqs2, alignment_t& aln);
double nuc_pi(Vector4d n, Vector5d pis);

#endif
