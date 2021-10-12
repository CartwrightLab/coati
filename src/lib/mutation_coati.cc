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
// const float omega = 0.2;  // github.com/reedacartwright/toycoati

namespace coati {
/**
 * \brief Create Muse \& Gaut (1994) substitution matrix.
 *
 * Given a branch length, create a 64x64 codon substitution P matrix based on
 *  Muse \& Gaut model. Using nucleotide substitution rates from Yang (1994).
 *
 * @param[in] br_len float branch length.
 * @param[in] omega float nonsynonymous-synonymous bias.
 * @param[in] nuc_freqs std::vector<coati::float_t> nucleotide frequencies
 *  (A,C,G,T).
 *
 * \return substitution P matrix (coati::Matrixf).
 */
coati::Matrixf mg94_p(float br_len, float omega,
                      const std::vector<coati::float_t>& nuc_freqs) {
    if(br_len <= 0) {
        throw std::out_of_range("Branch length must be positive.");
    }

    // Yang (1994) estimating the pattern of nucleotide substitution
    float nuc_q[4][4] = {{-0.818, 0.132, 0.586, 0.1},
                         {0.221, -1.349, 0.231, 0.897},
                         {0.909, 0.215, -1.322, 0.198},
                         {0.1, 0.537, 0.128, -0.765}};

    // MG94 model - doi:10.1534/genetics.108.092254
    Matrix64f Q = Matrix64f::Zero();
    float Pi[64];
    float w{NAN}, d = 0.0f;
    int x = 0, y = 0;
    uint8_t first_cod_mask = 48, second_cod_mask = 12, third_cod_mask = 3;

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

        Pi[i] = nuc_freqs[((i & 48) >> 4)] * nuc_freqs[((i & 12) >> 2)] *
                nuc_freqs[(i & 3)];
        float rowSum = 0.0;
        for(uint8_t j = 0; j < 64; j++) {
            if(i == j) {
                Q(i, j) = 0;
            } else if(coati::utils::cod_distance(i, j) > 1) {
                Q(i, j) = 0;
            } else {
                w = ((amino_group_table[i] == amino_group_table[j]) ? 1
                                                                    : omega);

                // split into cases to avoid use of pow (speed-up)
                if((i & first_cod_mask) != (j & first_cod_mask)) {
                    x = (i & first_cod_mask) >> 4;
                    y = (j & first_cod_mask) >> 4;
                } else if((i & second_cod_mask) != (j & second_cod_mask)) {
                    x = (i & second_cod_mask) >> 2;
                    y = (j & second_cod_mask) >> 2;
                } else if((i & third_cod_mask) != (j & third_cod_mask)) {
                    x = i & third_cod_mask;
                    y = j & third_cod_mask;
                }

                Q(i, j) = w * nuc_q[x][y];
            }
            rowSum += Q(i, j);
        }
        Q(i, i) = -rowSum;
        d += Pi[i] * rowSum;
    }

    // normalize
    Q = Q / d;

    Q = Q * br_len;
    Q = Q.exp();

    coati::Matrixf P(64, 64, Q);

    return P;
}

/// @private
TEST_CASE("mg94_p") {
    coati::Matrixf P(mg94_p(0.0133, 0.2, {0.308, 0.185, 0.199, 0.308}));

    for(int i = 0; i < 64; i++) {
        for(int j = 0; j < 64; j++) {
            CHECK(P(i, j) == doctest::Approx(mg94P[i][j]));
        }
    }
}

/**
 * \brief Create marginal 192x4 substitution P matrix give a 64x64
 *  substitution matrix.
 *
 * @param[in] P coati::Matrixf 64x64 codon substitution matrix.
 * @param[in] pi std::vector<coati::float_t> nucleotide frequencies (A,C,G,T).
 *
 * \return marginal 192x4 substitution matrix (coati::Matrixf).
 */
coati::Matrixf marginal_p(const coati::Matrixf& P,
                          const std::vector<coati::float_t>& pi) {
    float marg{NAN};

    coati::Matrixf p(192, 4);

    for(int cod = 0; cod < 64; cod++) {
        for(int pos = 0; pos < 3; pos++) {
            for(int nuc = 0; nuc < 4; nuc++) {
                marg = 0.0;
                switch(pos) {  // divide cases into each value of pos for speed
                               // up (reduce use of pow())
                case 0:
                    for(uint8_t i = 0; i < 64; i++) {
                        marg += (((i & 48) >> 4) == nuc ? P(cod, i) : 0.0f);
                    }
                    break;
                case 1:
                    for(uint8_t i = 0; i < 64; i++) {
                        marg += (((i & 12) >> 2) == nuc ? P(cod, i) : 0.0f);
                    }
                    break;
                case 2:
                    for(uint8_t i = 0; i < 64; i++) {
                        marg += ((i & 3) == nuc ? P(cod, i) : 0.0f);
                    }
                    break;
                }
                p(cod * 3 + pos, nuc) = ::logf(marg / pi[nuc]);
            }
        }
    }

    return p;
}

/// @private
TEST_CASE("marginal_p") {
    std::vector<coati::float_t> pi{0.308, 0.185, 0.199, 0.308};
    coati::Matrixf P = mg94_p(0.0133, 0.2, pi);
    coati::Matrixf p_marg = marginal_p(P, pi);

    for(int cod = 0; cod < 64; cod++) {
        for(int pos = 0; pos < 3; pos++) {
            float val = 0.f;
            for(int nuc = 0; nuc < 4; nuc++) {
                val += ::expf(p_marg(cod * 3 + pos, nuc)) * pi[nuc];
            }
            CHECK(val == doctest::Approx(1));  // sum per pos (all nuc) is 1
        }
    }
}
}  // namespace coati
