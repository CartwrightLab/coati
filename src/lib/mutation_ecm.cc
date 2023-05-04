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

#include <doctest/doctest.h>

#include <coati/mutation_ecm.hpp>

namespace coati {

/**
 * @brief Calculate number of transitions and transversion between two codons.
 *
 * @details Calculate number of ts and tv between two codons one position at a
 * time. Nucleotides are encoded as unsigned ints - A=0, C=1, G=2, T=3. If both
 * are odd (C,T) then it's at transversion, otherwhise (even - A,G) it's a
 * transition.
 *
 * @param[in] c1 uint8_t encoded codon (AAA=0,...,TTT=63).
 * @param[in] c2 uint8_t encoded codon.
 * @param[in,out] nts int number of transitions.
 * @param[in,out] ntv int number of transversions.
 *
 */
void nts_ntv(uint8_t c1, uint8_t c2, int& nts, int& ntv) {
    using coati::utils::get_nuc;
    nts = ntv = 0;
    if(c1 == c2) {
        return;  // if same codon, 0 ts & 0 tv - return
    }

    for(int i = 0; i < 3; i++) {
        uint8_t nt1 = coati::utils::get_nuc(c1, i);
        uint8_t nt2 = coati::utils::get_nuc(c2, i);

        if(nt1 == nt2) {
            continue;  // if same nuc, 0 ts & 0 tv - go to next nuc
        }
        ((nt1 % 2 == nt2 % 2) ? nts : ntv) += 1;
    }
}

/// @private
// GCOVR_EXCL_START
TEST_CASE("nts_ntv") {
    int nts = 0, ntv = 0;

    nts_ntv(0, 0, nts, ntv);  // AAA -> AAA
    CHECK_EQ(nts, 0);
    CHECK_EQ(ntv, 0);

    nts_ntv(0, 1, nts, ntv);  // AAA -> AAC
    CHECK_EQ(nts, 0);
    CHECK_EQ(ntv, 1);

    nts_ntv(39, 57, nts, ntv);  // GCT -> TTA
    CHECK_EQ(nts, 1);
    CHECK_EQ(ntv, 2);

    nts_ntv(21, 42, nts, ntv);  // CCC -> GGG
    CHECK_EQ(nts, 0);
    CHECK_EQ(ntv, 3);

    nts_ntv(42, 0, nts, ntv);  // GGG -> AAA
    CHECK_EQ(nts, 3);
    CHECK_EQ(ntv, 0);
}
// GCOVR_EXCL_STOP

/**
 * @brief Transition-transverison bias function.
 *
 * @details Used in the empirical codon model (ECM), k represents the the
 * relative strength of the transitions-transversion bias, modeled as a function
 * that depends on the number of transitions and transversion between two
 * codons. More information here https://doi.org/10.1093/molbev/msm064.
 *
 * @param[in] c1 uint8_t encoded codon (AAA=0, ... , TTT=63)
 * @param[in] c2 uint8_t encoded codon.
 * @param[in] model int model assumed to best fit the data.
 * @param[in] kappa float measure relative to the value implicit in
 *  exchangeabilities.
 *
 * @retval float relative strength of transition-transversion bias.
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

/// @private
// GCOVR_EXCL_START
TEST_CASE("k") {
    CHECK_EQ(k(0, 0, 0), 1);         // AAA -> AAA, ECM+f+omega
    CHECK_EQ(k(32, 0, 0), 1);        // GAA -> CTC, ECM+f+omega
    CHECK_EQ(k(47, 38), 1);          // GTT -> GCT, ECM+f+omega
    CHECK_EQ(k(22, 19), 1);          // CCG -> CAT, ECM+f+omega
    CHECK_EQ(k(0, 42, 1), 15.625);   // AAA -> GGG, ECM+F+omega+1k(ts)
    CHECK_EQ(k(32, 29, 1), 1);       // GAA -> CTC, ECM+F+omega+1k(ts)
    CHECK_EQ(k(47, 38, 1), 2.5);     // GTT -> GCT, ECM+F+omega+1k(ts)
    CHECK_EQ(k(21, 49, 1), 6.25);    // CCC -> TAT, ECM+F+omega+1k(ts)
    CHECK_EQ(k(0, 0, 2), 1);         // AAA -> AAA, ECM+F+omega+1k(tv)
    CHECK_EQ(k(32, 29, 2), 15.625);  // GAA -> CTC, ECM+F+omega+1k(tv)
    CHECK_EQ(k(47, 38, 2), 2.5);     // GTT -> GCT, ECM+F+omega+1k(tv)
    CHECK_EQ(k(22, 19, 2), 6.25);    // CCG -> CAT, ECM+F+omega+1k(tv)
}
// GCOVR_EXCL_STOP

/**
 * @brief Create Empirical Codon Model substitution P matrix.
 *
 * @param[in] br_len float brach length.
 * @param[in] omega float nonsynonymous-synonymous bias.
 *
 * @retval coati::Matrixf empirical codon model substitution P matrix.
 */
coati::Matrixf ecm_p(float br_len, float omega) {
    if(br_len <= 0) {
        throw std::out_of_range("Branch length must be positive.");
    }

    Matrix61f Q = Matrix61f::Zero();

    float d = 0.0;

    for(uint8_t i = 0; i < 61; i++) {
        float rowSum = 0.0;
        for(uint8_t j = 0; j < 61; j++) {
            // if codons i or j are stop codons value is zero (default)
            if(i == j) {
                continue;
            }
            if(amino_group[i] == amino_group[j]) {
                // same aminoacid group
                Q(i, j) = exchang[i][j] * ecm_pi[j] * k(i, j, 0);
            } else {
                // different aminoacid group
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

    return {61, 61, Q};
}

/**
 * @brief Create Empirical Codon Model (Kosiol et al. 2007) FST.
 *
 * @param[in] br_len float branch length;
 * @param[in] omega float nonsynonymous-synonymous bias.
 *
 * @retval coati::VectorFstStdArc empirical codon model FST.
 */
VectorFstStdArc ecm(float br_len, float omega) {
    coati::Matrixf P = ecm_p(br_len, omega);
    using coati::utils::get_nuc;

    // Add state 0 and make it the start state
    VectorFstStdArc ecm;
    ecm.AddState();
    ecm.SetStart(0);

    // Create FST
    int r = 1;
    for(uint8_t i = 0; i < 61; i++) {
        for(uint8_t j = 0; j < 61; j++) {
            add_arc(ecm, 0, r, get_nuc(i, 0) + 1, get_nuc(j, 0) + 1, P(i, j));
            add_arc(ecm, r, r + 1, get_nuc(i, 1) + 1, get_nuc(j, 1) + 1);
            add_arc(ecm, r + 1, 0, get_nuc(i, 2) + 1, get_nuc(j, 2) + 1);
            r = r + 2;
        }
    }

    // Set final state
    ecm.SetFinal(0, 0.0);

    return optimize(ecm);
}
}  // namespace coati
