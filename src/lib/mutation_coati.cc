/*
# Copyright (c) 2020-2021 Juan J. Garcia Mesa <juanjosegarciamesa@gmail.com>
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

#include <doctest/doctest.h>

#include <coati/mutation_coati.hpp>

/* nonsynonymous-synonymous bias (\omega) */
const float omega = 0.2;  // github.com/reedacartwright/toycoati

/* Muse & Gaut Model (1994) P matrix */
void mg94_q(Matrix64f& Q) {
    // Yang (1994) estimating the pattern of nucleotide substitution
    float nuc_freqs[4] = {0.308, 0.185, 0.199, 0.308};
    float nuc_q[4][4] = {{-0.818, 0.132, 0.586, 0.1},
                         {0.221, -1.349, 0.231, 0.897},
                         {0.909, 0.215, -1.322, 0.198},
                         {0.1, 0.537, 0.128, -0.765}};

    // MG94 model - doi:10.1534/genetics.108.092254
    Q = Matrix64f::Zero();
    Matrix64f Pi;
    float w{NAN}, d = 0.0f;
    int x = 0, y = 0;

    // construct transition matrix
    for(uint8_t i = 0; i < 64; i++) {
        // uint8_t codon;
        // (codon & 48) >> 4 = nt4_table encoding of first codon nucleotide
        // (codon & 12) >> 2 = nt4_table encoding of second codon nucleotide
        // (codon & 03) = nt4_table encoding of third codon nucleotide
        // e.g. 00 00 11 10 = 00 A T G = codon "ATG"
        // (00001110 & 48) >> 4 = (00001110 & 00110000) >> 4 = 00000000 >> 4 = 0
        // (A) (00001110 & 12) >> 2 = (00001110 & 00001100) >> 2 = 00001100 >> 2
        // = 3 (T)
        // (00001110 & 03) 		= (00001110 & 00000011) 	 =
        // 00000010 = 2 (G)

        Pi(i) = nuc_freqs[((i & 48) >> 4)] * nuc_freqs[((i & 12) >> 2)] *
                nuc_freqs[(i & 3)];
        float rowSum = 0.0;
        for(uint8_t j = 0; j < 64; j++) {
            if(i == j) {
                Q(i, j) = 0;
            } else if(cod_distance(i, j) > 1) {
                Q(i, j) = 0;
            } else {
                w = ((nt4_table[i] == nt4_table[j]) ? 1 : omega);

                // split into cases to avoid use of pow (speed-up)
                if((i & static_cast<uint8_t>(48)) !=
                   (j & static_cast<uint8_t>(48))) {
                    x = (i & static_cast<uint8_t>(48)) >> 4;
                    y = (j & static_cast<uint8_t>(48)) >> 4;
                } else if((i & static_cast<uint8_t>(12)) !=
                          (j & static_cast<uint8_t>(12))) {
                    x = (i & static_cast<uint8_t>(12)) >> 2;
                    y = (j & static_cast<uint8_t>(12)) >> 2;
                } else if((i & static_cast<uint8_t>(3)) !=
                          (j & static_cast<uint8_t>(3))) {
                    x = i & static_cast<uint8_t>(3);
                    y = j & static_cast<uint8_t>(3);
                }

                Q(i, j) = w * nuc_q[x][y];
            }
            rowSum += Q(i, j);
        }
        Q(i, i) = -rowSum;
        d += Pi(i) * rowSum;
    }

    // normalize
    Q = Q / d;
}

TEST_CASE("[mutation_coati.cc] mg94_q") {
    Matrix64f Q;
    mg94_q(Q);

    for(int i = 0; i < 64; i++) {
        CHECK(Q(i, i) < 0);  // main diagonal has negative values
        CHECK(Q.row(i).sum() == doctest::Approx(0));  // sum of every row is 0
    }

    CHECK((Q.array() == 0.0).count() ==
          64 * 64 - 64 * 10);  // only 10 values per row are non-zero

    // convert Q matrix to column vector
    Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Q2(Q);
    Eigen::Map<Eigen::VectorXf> mg94Q_vector(Q2.data(), 4096, 1);

    // check non-zero values
    for(int i = 0; i < 640; i++) {
        CHECK(mg94Q_vector[indexes[i]] == doctest::Approx(mg94Q[i]));
    }
}

/* Muse & Gaut Model (1994) P matrix given rate matrix and branch lenght */
void mg94_p(Matrix64f& P, float brlen) {
    if(brlen <= 0) {
        throw std::out_of_range("Branch length must be positive.");
    }

    Matrix64f Q;
    mg94_q(Q);

    Q = Q * brlen;
    P = Q.exp();
}

TEST_CASE("[mutation_coati.cc] mg94_p") {
    Matrix64f P;
    float branch_length = 0.0133;
    mg94_p(P, branch_length);

    for(int i = 0; i < 64; i++) {
        CHECK(P.row(i).sum() == doctest::Approx(1));  // sum of every row is 1
    }

    CHECK(!(P.array() < 0).any());  // all values are positive

    for(int i = 0; i < 64; i++) {
        for(int j = 0; j < 64; j++) {
            CHECK(P(i, j) == doctest::Approx(mg94P[i][j]));
        }
    }
}

/* Create marginal Muse and Gaut codon model P matrix*/
void mg94_marginal_p(Eigen::Tensor<float, 3>& p, Matrix64f& P) {
    float marg{NAN};

    for(int cod = 0; cod < 64; cod++) {
        for(int pos = 0; pos < 3; pos++) {
            for(int nuc = 0; nuc < 4; nuc++) {
                marg = 0.0;
                switch(pos) {  // divide cases into each value of pos for speed
                               // up (reduce use of pow())
                case 0:
                    for(uint8_t i = 0; i < 64; i++) {
                        marg += (((i & static_cast<uint8_t>(48)) >> 4) == nuc
                                     ? P(cod, i)
                                     : 0.0f);
                    }
                    break;
                case 1:
                    for(uint8_t i = 0; i < 64; i++) {
                        marg += (((i & static_cast<uint8_t>(12)) >> 2) == nuc
                                     ? P(cod, i)
                                     : 0.0f);
                    }
                    break;
                case 2:
                    for(uint8_t i = 0; i < 64; i++) {
                        marg +=
                            ((i & static_cast<uint8_t>(3)) == nuc ? P(cod, i)
                                                                  : 0.0f);
                    }
                    break;
                }
                p(cod, pos, nuc) = marg;
            }
        }
    }
}

TEST_CASE("[mutation_coati.cc] mg94_marginal_p") {
    Eigen::Tensor<float, 3> p(64, 3, 4);
    Matrix64f P;
    float branch_length = 0.0133;
    mg94_p(P, branch_length);
    mg94_marginal_p(p, P);

    for(int cod = 0; cod < 64; cod++) {
        for(int pos = 0; pos < 3; pos++) {
            float val = 0;
            for(int nuc = 0; nuc < 4; nuc++) {
                val += p(cod, pos, nuc);
                CHECK(!(p(cod, pos, nuc) < 0));
            }
            CHECK(val ==
                  doctest::Approx(1));  // sum per each pos (all nuc) is 1
        }
    }
}
