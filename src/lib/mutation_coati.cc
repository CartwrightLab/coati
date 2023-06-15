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

#include <doctest/doctest.h>

#include <coati/mutation_coati.hpp>

namespace coati {
/**
 * @brief Create Muse \& Gaut (1994) substitution matrix.
 *
 * @details Given a branch length, create a 61x61 codon substitution P matrix
 * based on Muse \& Gaut model. Using either nucleotide substitution rates from
 * Yang (1994) or the General Time Reversible (GTR) model from Tavaré (1986).
 *
 * @param[in] br_len float branch length.
 * @param[in] omega float nonsynonymous-synonymous bias.
 * @param[in] nuc_freqs std::vector<coati::float_t> nucleotide frequencies
 *  (A,C,G,T).
 * @param[in] sigma std::vector<coati::float_t> transition probabilities for GTR
 * substitution model.
 *
 * @retval coati::Matrixf substitution P matrix.
 */
coati::Matrixf mg94_p(float br_len, float omega,
                      const std::vector<coati::float_t>& nuc_freqs,
                      const std::vector<coati::float_t>& sigma) {
    using coati::utils::get_nuc;
    if(br_len <= 0) {
        throw std::out_of_range("Branch length must be positive.");
    }

    coati::Matrixf nuc_q(4, 4);

    if(std::any_of(sigma.cbegin(), sigma.cend(),
                   [](coati::float_t f) { return f > 0.f; })) {
        // Use GTR model for nuc_q
        nuc_q = gtr_q(nuc_freqs, sigma);
    } else {
        // Use Yang (1994) estimating the pattern of nucleotide substitution
        nuc_q = {{-0.818, 0.132, 0.586, 0.1},
                 {0.221, -1.349, 0.231, 0.897},
                 {0.909, 0.215, -1.322, 0.198},
                 {0.1, 0.537, 0.128, -0.765}};
    }

    // MG94 model - doi:10.1534/genetics.108.092254
    Matrix61f Q = Matrix61f::Zero();
    float Pi[61];
    float w{NAN}, d = 0.0f;
    int x = 0, y = 0;

    // construct transition matrix
    for(uint8_t i = 0; i < 61; i++) {
        Pi[i] = nuc_freqs[get_nuc(i, 0)] * nuc_freqs[get_nuc(i, 1)] *
                nuc_freqs[get_nuc(i, 2)];
        /* Codon frequency by multiplying nucleotide frequencies.
         *  Accessing nucleotide frequencies using the following:
         * (codon & 48) >> 4 = nt16_table encoding of first codon nucleotide
         * (codon & 12) >> 2 = nt16_table encoding of second codon nucleotide
         * (codon & 03) = nt16_table encoding of third codon nucleotide
         * e.g. 00 11 10 = A T G = codon "ATG" (note: omitting first 2 bits)
         * (001110 & 48) >> 4 =(001110 & 00110000) >> 4 = 000000 >> 4 = 0 (A)
         * (001110 & 12) >> 2 = (001110 & 001100) >> 2  = 001100 >> 2 = 3 (T)
         * (001110 & 03)      = (001110 & 000011)       = 00000010    = 2 (G)
         */

        float rowSum = 0.0;
        for(uint8_t j = 0; j < 61; j++) {
            if(i == j) {
                Q(i, j) = 0;
            } else if(coati::utils::cod_distance(i, j) > 1) {
                Q(i, j) = 0;
            } else {
                w = ((amino_group[i] == amino_group[j]) ? 1 : omega);

                // x,y = positions of two nucleotides involved in substitution
                if(get_nuc(i, 0) != get_nuc(j, 0)) {
                    x = get_nuc(i, 0);
                    y = get_nuc(j, 0);
                } else if(get_nuc(i, 1) != get_nuc(j, 1)) {
                    x = get_nuc(i, 1);
                    y = get_nuc(j, 1);
                } else {
                    x = get_nuc(i, 2);
                    y = get_nuc(j, 2);
                }

                Q(i, j) = w * nuc_q(x, y);
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

    return {61, 61, Q};
}

/// @private
// GCOVR_EXCL_START
TEST_CASE("mg94_p") {
    SUBCASE("default") {
        coati::Matrixf P(mg94_p(0.0133, 0.2, {0.308, 0.185, 0.199, 0.308}));

        for(int i = 0; i < 61; i++) {
            for(int j = 0; j < 61; j++) {
                CHECK_EQ(P(i, j), doctest::Approx(mg94P[i][j]));
            }
        }
    }
    SUBCASE("invalid branch length - fail") {
        CHECK_THROWS_AS(mg94_p(0.f, 0.2, {0.308, 0.185, 0.199, 0.308}),
                        std::out_of_range);
        CHECK_THROWS_AS(mg94_p(-0.02, 0.2, {0.308, 0.185, 0.199, 0.308}),
                        std::out_of_range);
    }
}
// GCOVR_EXCL_STOP

/**
 * @brief Create marginal 183x15 substitution P matrix.
 *
 * @details Give a 61x61 codon substitution matrix, marginalized by codon,
 * phase, and nucleotide (i.e. P(nuc | codon, phase), where phase = {0,1,2}).
 * Dimensions are 61 codons * 3 phases * 15 nucleotides = 183x15.
 *  note: 15 nucleotides defined in IUPAC code
 *        https://www.bioinformatics.org/sms/iupac.html.
 *
 * @param[in] P coati::Matrixf 61x61 codon substitution matrix.
 * @param[in] pi std::vector<coati::float_t> nucleotide frequencies (A,C,G,T).
 * @param[in] amb coati::AmbiguousNucs how to calculate probabilities for
 * ambiguous nucleotides.
 *
 * @retval coati::Matrixf marginal 183x4 substitution matrix.
 */
coati::Matrixf marginal_p(const coati::Matrixf& P,
                          const std::vector<coati::float_t>& pi,
                          const coati::AmbiguousNucs amb) {
    using coati::utils::get_nuc;
    float marg{NAN};

    coati::Matrixf p(183, 15);

    for(size_t cod = 0; cod < P.rows(); cod++) {
        for(int nuc = 0; nuc < 4; nuc++) {  // A,C,G,T
            // position 0
            marg = 0.0;
            for(uint8_t i = 0; i < P.cols(); i++) {
                marg += (get_nuc(i, 0) == nuc ? P(cod, i) : 0.0f);
            }
            p(cod * 3, nuc) = ::logf(marg / pi[nuc]);

            // position 1
            marg = 0.0;
            for(uint8_t i = 0; i < P.cols(); i++) {
                marg += (get_nuc(i, 1) == nuc ? P(cod, i) : 0.0f);
            }
            p(cod * 3 + 1, nuc) = ::logf(marg / pi[nuc]);

            // position 2
            marg = 0.0;
            for(uint8_t i = 0; i < P.cols(); i++) {
                marg += (get_nuc(i, 2) == nuc ? P(cod, i) : 0.0f);
            }
            p(cod * 3 + 2, nuc) = ::logf(marg / pi[nuc]);
        }
    }

    if(amb == AmbiguousNucs::AVG) {
        ambiguous_sum_p(p);
    } else {
        ambiguous_best_p(p);
    }

    return p;
}

/// @private
// GCOVR_EXCL_START
TEST_CASE("marginal_p") {
    std::vector<coati::float_t> pi{0.308, 0.185, 0.199, 0.308};
    coati::Matrixf P = mg94_p(0.0133, 0.2, pi);
    coati::Matrixf p_marg = marginal_p(P, pi, AmbiguousNucs::AVG);

    for(int cod = 0; cod < 61; cod++) {
        for(int pos = 0; pos < 3; pos++) {
            // NOLINTNEXTLINE(clang-diagnostic-unused-but-set-variable)
            float val = 0.f;
            for(int nuc = 0; nuc < 4; nuc++) {
                val += ::expf(p_marg(cod * 3 + pos, nuc)) * pi[nuc];
            }
            CHECK_EQ(val, doctest::Approx(1));  // sum per pos (all nuc) is 1
        }
    }
}
// GCOVR_EXCL_STOP

/**
 * @brief Probabilities for ambiguous nucs in marginal model using averages.
 *
 * @details The probability for each nucleotide is calculated by averaging over
 * all possibilities (e.g. N = avg(A, C, G, T); R = avg(A, G)).
 *
 * @param[in,out] p coati::Matrixf marginal substitution matrix.
 *
 */
void ambiguous_sum_p(coati::Matrixf& p) {
    using coati::utils::log_sum_exp;
    size_t row{0};
    for(int cod = 0; cod < 61; cod++) {
        for(size_t pos = 0; pos < 3; pos++) {
            row = cod * 3 + pos;
            // R: purine       A or G
            p(row, 4) = log_sum_exp(p(row, 0), p(row, 2));
            // Y: pyrymidine   C or T
            p(row, 5) = log_sum_exp(p(row, 1), p(row, 3));
            // M: amino group  A or C
            p(row, 6) = log_sum_exp(p(row, 0), p(row, 1));
            // K: keto group   G or T
            p(row, 7) = log_sum_exp(p(row, 2), p(row, 3));
            // S: strong inter C or G
            p(row, 8) = log_sum_exp(p(row, 1), p(row, 2));
            // W: weak interac A or T
            p(row, 9) = log_sum_exp(p(row, 0), p(row, 3));
            // B: not A
            p(row, 10) =
                log_sum_exp(log_sum_exp(p(row, 1), p(row, 2)), p(row, 3));
            // D: not C
            p(row, 11) =
                log_sum_exp(log_sum_exp(p(row, 0), p(row, 2)), p(row, 3));
            // H: not G
            p(row, 12) =
                log_sum_exp(log_sum_exp(p(row, 0), p(row, 1)), p(row, 3));
            // V: not T
            p(row, 13) =
                log_sum_exp(log_sum_exp(p(row, 0), p(row, 1)), p(row, 2));
            // N: any
            p(row, 14) = log_sum_exp(
                log_sum_exp(log_sum_exp(p(row, 0), p(row, 1)), p(row, 2)),
                p(row, 3));
        }
    }
}

/**
 * @brief Probabilities for ambiguous nucs in marginal model taking best prob.
 *
 * @details The probability for each nucleotide is calculated by taking the best
 * (highest) of all possibilities (e.g. N = max(A, C, G, T); R = max(A, G)).
 *
 * @param[in,out] p coati::Matrixf marginal substitution matrix.
 *
 */
void ambiguous_best_p(coati::Matrixf& p) {
    size_t row{0};
    for(int cod = 0; cod < 61; cod++) {
        for(size_t pos = 0; pos < 3; pos++) {
            row = cod * 3 + pos;

            p(row, 4) =
                std::max(p(row, 0), p(row, 2));  // R: purine       A or G
            p(row, 5) =
                std::max(p(row, 1), p(row, 3));  // Y: pyrimidine   C or T
            p(row, 6) =
                std::max(p(row, 0), p(row, 1));  // M: amino group  A or C
            p(row, 7) =
                std::max(p(row, 2), p(row, 3));  // K: keto group   G or T
            p(row, 8) =
                std::max(p(row, 1), p(row, 2));  // S: strong inter C or G
            p(row, 9) =
                std::max(p(row, 0), p(row, 3));  // W: weak interac A or T
            p(row, 10) =
                std::max({p(row, 1), p(row, 2), p(row, 3)});  // B: not A
            p(row, 11) =
                std::max({p(row, 0), p(row, 2), p(row, 3)});  // D: not C
            p(row, 12) =
                std::max({p(row, 0), p(row, 1), p(row, 3)});  // H: not G
            p(row, 13) =
                std::max({p(row, 0), p(row, 1), p(row, 2)});  // V: not T
            p(row, 14) = std::max(
                {p(row, 0), p(row, 1), p(row, 2), p(row, 3)});  // N: any
        }
    }
}

/**
 * @brief Create GTR substitution model matrix.
 *
 * @param[in] pi std::vector<coati::float_t> nucleotide frequencies.
 * @param[in] sigma std::vector<coati::float_t> sigma parameters (6) for GTR
 * model.
 *
 * @retval coati::Matrixf GTR Q matrix.
 */
coati::Matrixf gtr_q(const std::vector<coati::float_t>& pi,
                     const std::vector<coati::float_t>& sigma) {
    //   |        A      |       C       |       G       |       T       |
    // A |        -      | pi_C*sigma_AC | pi_G*sigma_AG | pi_T*sigma_AT |
    // C | pi_A*sigma_AC |        -      | pi_G*sigma_CG | pi_T*sigma_CT |
    // G | pi_A*sigma_AG | pi_C*sigma_GC |       -       | pi_T*sigma_GT |
    // T | pi_A*sigma_AT | pi_C*sigma_CT | pi_G*sigma_GT |       -       |

    if(std::any_of(sigma.cbegin(), sigma.cend(),
                   [](coati::float_t f) { return f < 0.f || f > 1.f; })) {
        throw std::invalid_argument("Sigma values must be in range [0,1].");
    }

    coati::Matrixf gtr_mat(4, 4);

    // set sigmas
    gtr_mat(0, 1) = gtr_mat(1, 0) = sigma[0];  // sigma_AC
    gtr_mat(0, 2) = gtr_mat(2, 0) = sigma[1];  // sigma_AG
    gtr_mat(0, 3) = gtr_mat(3, 0) = sigma[2];  // sigma_AT
    gtr_mat(1, 2) = gtr_mat(2, 1) = sigma[3];  // sigma_GC
    gtr_mat(1, 3) = gtr_mat(3, 1) = sigma[4];  // sigma_CT
    gtr_mat(2, 3) = gtr_mat(3, 2) = sigma[5];  // sigma_GT

    // multiply by corresponding pi
    for(size_t i = 0; i < 4; i++) {
        for(size_t j = 0; j < 4; j++) {
            gtr_mat(i, j) *= pi[j];
        }
    }

    // set major diagonal
    gtr_mat(0, 0) = -(gtr_mat(0, 1) + gtr_mat(0, 2) + gtr_mat(0, 3));
    gtr_mat(1, 1) = -(gtr_mat(1, 0) + gtr_mat(1, 2) + gtr_mat(1, 3));
    gtr_mat(2, 2) = -(gtr_mat(2, 0) + gtr_mat(2, 1) + gtr_mat(2, 3));
    gtr_mat(3, 3) = -(gtr_mat(3, 0) + gtr_mat(3, 1) + gtr_mat(3, 2));

    return gtr_mat;
}

/// @private
// GCOVR_EXCL_START
TEST_CASE("gtr_q") {
    coati::Matrixf gtr(gtr_q({0.308, 0.185, 0.199, 0.308},
                             {0.009489730, 0.039164824, 0.004318182,
                              0.015438693, 0.038734091, 0.008550000}));

    coati::Matrixf gtr_correct = {
        {-0.010879400, 0.001755600, 0.00779380, 0.00133000},
        {0.002922837, -0.017925237, 0.00307230, 0.01193010},
        {0.012062766, 0.002856158, -0.01755232, 0.00263340},
        {0.001330000, 0.007165807, 0.00170145, -0.01019726}};

    for(int i = 0; i < 4; i++) {
        for(int j = 0; j < 4; j++) {
            CHECK_EQ(gtr(i, j), doctest::Approx(gtr_correct(i, j)));
        }
    }

    SUBCASE("Sigma values out of range") {
        REQUIRE_THROWS_AS(gtr_q({0.308, 0.185, 0.199, 0.308},
                                {-0.009489730, 0.039164824, 0.004318182,
                                 0.015438693, 0.038734091, 0.008550000}),
                          std::invalid_argument);

        REQUIRE_THROWS_AS(gtr_q({0.308, 0.185, 0.199, 0.308},
                                {0.009489730, 0.039164824, 0.004318182,
                                 0.015438693, 1.038734091, 0.008550000}),
                          std::invalid_argument);
    }
}
// GCOVR_EXCL_STOP
}  // namespace coati
