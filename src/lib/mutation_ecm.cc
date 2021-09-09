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

#include <doctest/doctest.h>

#include <coati/mutation_ecm.hpp>

/* nonsynonymous-synonymous bias (\omega) */
// const float omega = 0.2;  // github.com/reedacartwright/toycoati
// const float kappa = 2.5;  // Kosiol et al. 2007, supplemental material

/* calculate number of transitions and transversions between codons c1 and c2*/
void nts_ntv(uint8_t c1, uint8_t c2, int& nts, int& ntv) {
    nts = ntv = 0;
    if(c1 == c2) return;
    for(int i = 0; i < 3; i++) {
        uint8_t nt1 =
            ((c1 & static_cast<uint8_t>(48 / pow(4, i))) >> (4 - 2 * i));
        uint8_t nt2 =
            ((c2 & static_cast<uint8_t>(48 / pow(4, i))) >> (4 - 2 * i));

        if(nt1 == nt2) continue;
        ((nt1 % 2 == nt2 % 2) ? nts : ntv) += 1;
    }
}

TEST_CASE("nts_ntv") {
    int nts = 0, ntv = 0;

    nts_ntv(0, 0, nts, ntv);  // AAA -> AAA
    CHECK(nts == 0);
    CHECK(ntv == 0);

    nts_ntv(0, 1, nts, ntv);  // AAA -> AAC
    CHECK(nts == 0);
    CHECK(ntv == 1);

    nts_ntv(39, 60, nts, ntv);  // GCT -> TTA
    CHECK(nts == 1);
    CHECK(ntv == 2);

    nts_ntv(21, 42, nts, ntv);  // CCC -> GGG
    CHECK(nts == 0);
    CHECK(ntv == 3);

    nts_ntv(42, 0, nts, ntv);  // GGG -> AAA
    CHECK(nts == 3);
    CHECK(ntv == 0);
}

/* transition-transversion bias function, depending on # of ts and tv (Nts,Ntv)
 */
float k(uint8_t c1, uint8_t c2, int model, float kappa) {
    int nts = 0, ntv = 0;
    nts_ntv(c1, c2, nts, ntv);
    switch(model) {
    case 0:
        return 1;  // ECM+f+omega. Assumes ts-tv bias is accounted for
    case 1:
        return powf(kappa, static_cast<float>(nts));  // ECM+F+omega+1k(ts)
    case 2:
        return powf(kappa, static_cast<float>(ntv));  // ECM+F+omega+1k(tv)
    default:
        return 1;  // ECM+f+omega. Assumes ts-tv bias is accounted for
    }

    return 0;
}

TEST_CASE("k") {
    CHECK(k(0, 0, 0) == 1);         // AAA -> AAA, ECM+f+omega
    CHECK(k(32, 0, 0) == 1);        // GAA -> CTC, ECM+f+omega
    CHECK(k(47, 38) == 1);          // GTT -> GCT, ECM+f+omega
    CHECK(k(22, 19) == 1);          // CCG -> CAT, ECM+f+omega
    CHECK(k(0, 42, 1) == 15.625);   // AAA -> GGG, ECM+F+omega+1k(ts)
    CHECK(k(32, 29, 1) == 1);       // GAA -> CTC, ECM+F+omega+1k(ts)
    CHECK(k(47, 38, 1) == 2.5);     // GTT -> GCT, ECM+F+omega+1k(ts)
    CHECK(k(21, 51, 1) == 6.25);    // CCC -> TAT, ECM+F+omega+1k(ts)
    CHECK(k(0, 0, 2) == 1);         // AAA -> AAA, ECM+F+omega+1k(tv)
    CHECK(k(32, 29, 2) == 15.625);  // GAA -> CTC, ECM+F+omega+1k(tv)
    CHECK(k(47, 38, 2) == 2.5);     // GTT -> GCT, ECM+F+omega+1k(tv)
    CHECK(k(22, 19, 2) == 6.25);    // CCG -> CAT, ECM+F+omega+1k(tv)
}

/* Empirical Codon Model P matrix */
coati::Matrixf ecm_p(float br_len, float omega) {
    if(br_len <= 0) {
        throw std::out_of_range("Branch length must be positive.");
    }

    Matrix64f Q = Matrix64f::Zero();

    float d = 0.0;

    for(uint8_t i = 0; i < 64; i++) {
        float rowSum = 0.0;
        for(uint8_t j = 0; j < 64; j++) {
            // check if codons i or j are stop codons
            if(i == j || amino_group_table[i] == '*' || amino_group_table[j] == '*') {
                continue;
            }
            if(amino_group_table[i] == amino_group_table[j]) {
                Q(i, j) = exchang[i][j] * ecm_pi[j] * k(i, j, 0);
            } else {  // amino_group_table[i] != amino_group_table[j]{
                Q(i, j) = exchang[i][j] * ecm_pi[j] * k(i, j, 0) * omega;
            }
            rowSum += Q(i, j);
        }
        Q(i, i) = -rowSum;
        d += ecm_pi[i] * rowSum;
    }

    // normalize
    Q = Q / d;

    // P matrix
    Q = Q * br_len;
    Q = Q.exp();

    coati::Matrixf P(64, 64, Q);

    return P;
}

/* Empirical Codon Model (Kosiol et al. 2007) FST */
VectorFstStdArc ecm(float br_len, float omega) {
    coati::Matrixf P = ecm_p(br_len, omega);

    // Add state 0 and make it the start state
    VectorFstStdArc ecm;
    ecm.AddState();
    ecm.SetStart(0);

    int r = 1;
    for(uint8_t i = 0; i < 64; i++) {
        for(uint8_t j = 0; j < 64; j++) {
            add_arc(ecm, 0, r, ((i & 48) >> 4) + 1, ((j & 48) >> 4) + 1,
                    P(i, j));
            add_arc(ecm, r, r + 1, ((i & 12) >> 2) + 1, ((j & 12) >> 2) + 1);
            add_arc(ecm, r + 1, 0, (i & 3) + 1, (j & 3) + 1);
            r = r + 2;
        }
    }

    // Set final state
    ecm.SetFinal(0, 0.0);
    return optimize(ecm);
}
