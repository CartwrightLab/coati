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

#include <coati/profile_aln.hpp>
#include <string>
#include <typeinfo>
#include <vector>

/* Create profile given a sequence */
Eigen::MatrixXf create_profile(std::string seq) {
    std::vector<std::string> vector_seq{std::move(seq)};
    return create_profile(vector_seq);
}

/* Create profile given an alignment */
Eigen::MatrixXf create_profile(std::vector<std::string>& aln) {
    for(const auto& s : aln) {
        if(s.length() != aln[0].length()) {
            throw std::invalid_argument(
                "Profile matrix requires all strings of same length.");
        }
    }

    auto cols = static_cast<int>(aln.at(0).length());
    auto rows = aln.size();
    Eigen::MatrixXf profile = Eigen::MatrixXf::Zero(4, cols);

    int nuc{0};
    for(int j = 0; j < cols; j++) {      // for each column
        for(int i = 0; i < rows; i++) {  // for each row
            switch(aln.at(i).at(j)) {
            case 'A':
            case 'a':
                nuc = 0;
                break;
            case 'C':
            case 'c':
                nuc = 1;
                break;
            case 'G':
            case 'g':
                nuc = 2;
                break;
            case 'T':
            case 't':
                nuc = 3;
                break;
            case '-':
                continue;
            }
            profile(nuc, j) += 1.0f / static_cast<float>(rows);
        }
    }

    return profile;
}

TEST_CASE("[utils.cc] create_profile") {
    SUBCASE("alignment") {
        std::vector<std::string> aln = {"CTCTGGATAGTG", "CT----ATAGTG",
                                        "CTCT---TAGTG", "CTCTG--TAGTG"};
        Eigen::MatrixXf profile = create_profile(aln);
        Eigen::MatrixXf result(4, 12);
        result << 0, 0, 0, 0, 0, 0, 0.5, 0, 1, 0, 0, 0, 1, 0, 0.75, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0.5, 0.25, 0, 0, 0, 1, 0, 1, 0, 1, 0,
            0.75, 0, 0, 0, 1, 0, 0, 1, 0;

        std::cout << profile << std::endl;
        CHECK(profile == result);
    }

    SUBCASE("sequence") {
        std::string seq = "CTCTGGATAGTG";
        Eigen::MatrixXf profile = create_profile(seq);
        Eigen::MatrixXf result(4, 12);
        result << 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1,
            0, 0, 1, 0;

        CHECK(profile == result);
    }
}

/* Return value for  */
float transition(Matrix4x3d cod, int pos, Vector4f nuc,
                 const Eigen::Tensor<float, 3>& p) {
    float val = 0;
    pos = pos == 0 ? 2 : --pos;

    // find highest value nucleotide per codon position
    Eigen::VectorXf::Index max_pos0 = 0, max_pos1 = 0, max_pos2 = 0;
    float max_cod0 = cod.col(0).maxCoeff(&max_pos0);
    float max_cod1 = cod.col(1).maxCoeff(&max_pos1);
    float max_cod2 = cod.col(2).maxCoeff(&max_pos2);

    // codon to int value (AAA->0, AAC->1, ... TTT-> 63)
    auto cod_index = (max_pos0 << 4) + (max_pos1 << 2) + max_pos2;

    // weighted average for each nuc freq with highest codon
    for(int i = 0; i < 4; i++) {
        val += max_cod0 * max_cod1 * max_cod2 * nuc[i] * p(cod_index, pos, i);
    }

    // get second highest nuc value per codon position
    cod.col(0)(max_pos0) = -1;
    cod.col(1)(max_pos1) = -1;
    cod.col(2)(max_pos2) = -1;
    Eigen::VectorXf::Index max2_pos0 = 0, max2_pos1 = 0, max2_pos2 = 0;
    float max2_cod0 = cod.col(0).maxCoeff(&max2_pos0);
    float max2_cod1 = cod.col(1).maxCoeff(&max2_pos1);
    float max2_cod2 = cod.col(2).maxCoeff(&max2_pos2);

    // codon: 2nd highest nuc, highest nuc, highest nuc
    cod_index = (max2_pos0 << 4) + (max_pos1 << 2) + max_pos2;
    for(int i = 0; i < 4; i++) {
        val += max2_cod0 * max_cod1 * max_cod2 * nuc[i] * p(cod_index, pos, i);
    }

    // codon: highest nuc, 2nd highest nuc, highest nuc
    cod_index = (max_pos0 << 4) + (max2_pos1 << 2) + max_pos2;
    for(int i = 0; i < 4; i++) {
        val += max_cod0 * max2_cod1 * max_cod2 * nuc[i] * p(cod_index, pos, i);
    }

    // codon: highest nuc, highest nuc, 2nd highest nuc
    cod_index = (max_pos0 << 4) + (max_pos1 << 2) + max2_pos2;
    for(int i = 0; i < 4; i++) {
        val += max_cod0 * max_cod1 * max2_cod2 * nuc[i] * p(cod_index, pos, i);
    }

    return val;
}

TEST_CASE("[profile_aln.cc] transition") {
    Vector4f nuc;
    Matrix4x3d cod;

    Matrix64f P;
    float br_len = 0.0133;
    mg94_p(P, br_len);
    Eigen::Tensor<float, 3> p(64, 3, 4);
    mg94_marginal_p(p, P);

    SUBCASE("codon aaa, all nucs") {
        cod << 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0;

        CHECK(transition(cod, 0, cod.col(0), p) == doctest::Approx(0.99352));
        CHECK(transition(cod, 1, cod.col(0), p) == doctest::Approx(0.99831));
        CHECK(transition(cod, 2, cod.col(0), p) == doctest::Approx(0.99832));
        nuc << 0, 1, 0, 0;
        CHECK(transition(cod, 0, nuc, p) == doctest::Approx(0.00027));
        CHECK(transition(cod, 1, nuc, p) == doctest::Approx(0.00027));
        CHECK(transition(cod, 2, nuc, p) == doctest::Approx(0.00027));
        nuc << 0, 0, 1, 0;
        CHECK(transition(cod, 0, nuc, p) == doctest::Approx(0.00599));
        CHECK(transition(cod, 1, nuc, p) == doctest::Approx(0.00121));
        CHECK(transition(cod, 2, nuc, p) == doctest::Approx(0.00121));
        nuc << 0, 0, 0, 1;
        CHECK(transition(cod, 0, nuc, p) == doctest::Approx(0.00021));
        CHECK(transition(cod, 1, nuc, p) == doctest::Approx(0.00021));
        CHECK(transition(cod, 2, nuc, p) == doctest::Approx(0.00021));
    }

    SUBCASE("codon ACA {0.4^3}, homogeneous nucleotide freq") {
        nuc << 0.25, 0.25, 0.25, 0.25;
        cod << 0.4, 0.2, 0.4, 0.3, 0.4, 0.1, 0.2, 0.3, 0.2, 0.1, 0.1, 0.3;

        CHECK(transition(cod, 0, nuc, p) == doctest::Approx(0.052));
        CHECK(transition(cod, 1, nuc, p) == doctest::Approx(0.052));
        CHECK(transition(cod, 2, nuc, p) == doctest::Approx(0.052));
    }

    SUBCASE("codon CGT {0.4,0.5,0.4}, nucleotide freq {0.3,0.2,0.2,0.3}") {
        nuc << 0.3, 0.2, 0.2, 0.3;
        cod << 0.3, 0.2, 0.3, 0.4, 0.2, 0.1, 0.2, 0.5, 0.2, 0.1, 0.1, 0.4;

        CHECK(transition(cod, 0, nuc, p) == doctest::Approx(0.06945));
        CHECK(transition(cod, 1, nuc, p) == doctest::Approx(0.05244));
        CHECK(transition(cod, 2, nuc, p) == doctest::Approx(0.04964));
    }
}

/* Gotoh dynamic programming alignment with marginal Muse & Gaut p matrix for
 * profile matrices  */
int gotoh_profile_marginal(std::vector<std::string> seqs1,
                           std::vector<std::string> seqs2, alignment_t& aln,
                           Matrix64f& P_m) {
    // P matrix for marginal Muse and Gaut codon model
    Eigen::Tensor<float, 3> p(64, 3, 4);
    mg94_marginal_p(p, P_m);

    Eigen::MatrixXf pro1 = create_profile(seqs1);
    Eigen::MatrixXf pro2 = create_profile(seqs2);

    int m = static_cast<int>(pro1.cols());
    int n = static_cast<int>(pro2.cols());

    // assert that length of 1st sequence (ref) is multiple of 3
    if(m % 3 != 0) {
        throw std::invalid_argument(
            "Reference CDS length must be of length multiple of 3.");
    }

    // DP matrices for match/mismatch (D), insertion (P), deletion (Q)
    Eigen::MatrixXf D = Eigen::MatrixXf::Constant(
        m + 1, n + 1, std::numeric_limits<float>::max());
    Eigen::MatrixXf P = Eigen::MatrixXf::Constant(
        m + 1, n + 1, std::numeric_limits<float>::max());
    Eigen::MatrixXf Q = Eigen::MatrixXf::Constant(
        m + 1, n + 1, std::numeric_limits<float>::max());

    // backtracking info matrices for match/mismatch (Bd), insert (Bp), and
    // deletion (Bq)
    Eigen::MatrixXi Bd = Eigen::MatrixXi::Constant(m + 1, n + 1, -1);
    Eigen::MatrixXi Bp = Eigen::MatrixXi::Constant(m + 1, n + 1, -1);
    Eigen::MatrixXi Bq = Eigen::MatrixXi::Constant(m + 1, n + 1, -1);

    float insertion = logf(0.001);
    float deletion = logf(0.001);
    float insertion_ext = logf(1.0f - (1.0 / 6.0));
    float deletion_ext = logf(1.0f - (1.0 / 6.0));
    float no_insertion = logf(1.0f - 0.001);
    float no_deletion = logf(1.0f - 0.001);
    float no_insertion_ext = logf(1.0f / 6.0);
    float no_deletion_ext = logf(1.0f / 6.0);

    Vector5f nuc_freqs;
    nuc_freqs << logf(0.308), logf(0.185), logf(0.199), logf(0.308), logf(0.25);

    // DP and backtracking matrices initialization

    // fill first values on D that are independent
    D(0, 0) = 0.0;  // 0.0;
    Bd(0, 0) = 0;
    D(0, 1) = -insertion - nuc_pi(pro2.col(0), nuc_freqs) - no_insertion_ext;
    P(0, 1) = -insertion - nuc_pi(pro2.col(0), nuc_freqs) - no_insertion_ext;
    Bd(0, 1) = 1;
    Bp(0, 1) = 2;
    D(1, 0) = -no_insertion - deletion - no_deletion_ext;
    Q(1, 0) = -no_insertion - deletion - no_deletion_ext;
    Bd(1, 0) = 2;
    Bq(1, 0) = 2;

    // fill first row of D
    if(n + 1 >= 2) {
        for(int j = 2; j < n + 1; j++) {
            D(0, j) = D(0, j - 1) - insertion_ext -
                      nuc_pi(pro2.col(j - 1), nuc_freqs);
            P(0, j) = P(0, j - 1) - insertion_ext -
                      nuc_pi(pro2.col(j - 1), nuc_freqs);
            Bd(0, j) = 1;
            Bp(0, j) = 1;
        }
    }

    // fill first column of D
    if(m + 1 >= 2) {
        for(int i = 2; i < m + 1; i++) {
            D(i, 0) = D(i - 1, 0) - deletion_ext;
            Q(i, 0) = Q(i - 1, 0) - deletion_ext;
            Bd(i, 0) = 2;
            Bq(i, 0) = 1;
        }
    }

    Matrix4x3d codon;
    float p1 = NAN, p2 = NAN, q1 = NAN, q2 = NAN, d = NAN;

    for(int i = 1; i < m + 1; i++) {
        // codon = seq_a.substr((((i-1)/3)*3),3); // current codon
        codon = pro1.block(0, (((i - 1) / 3) * 3), 4, 3);  // 3x3
        for(int j = 1; j < n + 1; j++) {
            // insertion
            p1 = P(i, j - 1) - insertion_ext -
                 nuc_pi(pro2.col(j - 1), nuc_freqs);
            p2 = Bd(i, j - 1) == 0
                     ? D(i, j - 1) - insertion -
                           nuc_pi(pro2.col(j - 1), nuc_freqs) - no_insertion_ext
                 : Bd(i, j - 1) == 1 ? D(i, j - 1) - insertion_ext -
                                           nuc_pi(pro2.col(j - 1), nuc_freqs)
                                     : std::numeric_limits<float>::max();
            P(i, j) = std::fmin(p1, p2);
            Bp(i, j) =
                p1 < p2
                    ? 1
                    : 2;  // 1 is insertion extension, 2 is insertion opening

            // deletion
            q1 = Q(i - 1, j) - deletion_ext;
            q2 = Bd(i - 1, j) == 0
                     ? D(i - 1, j) - no_insertion - deletion - no_deletion_ext
                 : Bd(i - 1, j) == 1 ? D(i - 1, j) - no_deletion_ext - deletion
                                     : D(i - 1, j) - deletion_ext;
            Q(i, j) = std::fmin(q1, q2);
            Bq(i, j) =
                q1 < q2 ? 1
                        : 2;  // 1 is deletion extension, 2 is deletion opening

            // match/mismatch
            if(Bd(i - 1, j - 1) == 0) {
                d = D(i - 1, j - 1) - no_insertion - no_deletion -
                    logf(transition(codon, (i) % 3, pro2.col(j - 1), p));
            } else if(Bd(i - 1, j - 1) == 1) {
                d = D(i - 1, j - 1) - no_deletion -
                    logf(transition(codon, (i) % 3, pro2.col(j - 1), p));
            } else {
                d = D(i - 1, j - 1) -
                    logf(transition(codon, (i) % 3, pro2.col(j - 1), p));
            }

            // D(i,j) = highest weight between insertion, deletion, and
            // match/mismatch
            //	in this case, lowest (-logf(weight)) value
            if(d < P(i, j)) {
                if(d < Q(i, j)) {
                    D(i, j) = d;
                    Bd(i, j) = 0;
                } else {
                    D(i, j) = Q(i, j);
                    Bd(i, j) = 2;
                }
            } else {
                if(P(i, j) < Q(i, j)) {
                    D(i, j) = P(i, j);
                    Bd(i, j) = 1;
                } else {
                    D(i, j) = Q(i, j);
                    Bd(i, j) = 2;
                }
            }
        }
    }

    aln.weight += D(m, n);  // weight

    // backtracking to obtain alignment
    return backtracking_profile(Bd, Bp, Bq, seqs1, seqs2, aln);
}

TEST_CASE("[profile_aln.cc] gotoh_profile_marginal") {
    std::vector<std::string> seqs1, seqs2;
    alignment_t aln, aln_pair;
    Matrix64f P;
    float branch = 0.0133;
    mg94_p(P, branch);

    SUBCASE("pairwise") {  // ensure result is identical to non-profile aln
        seqs1 = {"CTCTGG"};
        seqs2 = {"CCTGG"};

        REQUIRE(gotoh_profile_marginal(seqs1, seqs2, aln, P) == 0);
        CHECK(aln.weight == doctest::Approx(8.73227));
        CHECK(aln.f.seq_data[0].compare("CTCTGG") == 0);
        CHECK(aln.f.seq_data[1].compare("C-CTGG") == 0);
        REQUIRE(mg94_marginal({"CTCTGG", "CCTGG"}, aln_pair, P) == 0);
        CHECK(aln.f.seq_data == aln_pair.f.seq_data);
        CHECK(aln.weight == aln_pair.weight);
    }
}

/* Backtrack dynamic programming alignment of profile matrices and retrieve aln
 */
int backtracking_profile(Eigen::MatrixXi Bd, Eigen::MatrixXi Bp,
                         Eigen::MatrixXi Bq, std::vector<std::string> seqs1,
                         std::vector<std::string> seqs2, alignment_t& aln) {
    int i = static_cast<int>(seqs1[0].length());
    int j = static_cast<int>(seqs2[0].length());
    std::vector<std::string> aln_seqs(seqs1.size() + seqs2.size());

    while((i != 0) || (j != 0)) {
        switch(Bd(i, j)) {
        case 0:  // match/mismatch
            for(size_t s = 0; s < seqs1.size(); s++) {
                aln_seqs[s].insert(0, 1, seqs1[s][i - 1]);
            }
            for(size_t s = 0; s < seqs2.size(); s++) {
                aln_seqs[seqs1.size() + s].insert(0, 1, seqs2[s][j - 1]);
            }
            i--;
            j--;
            break;
        case 1:  // insertion
            while(Bp(i, j) == 1) {
                for(size_t s = 0; s < seqs1.size(); s++) {
                    aln_seqs[s].insert(0, 1, '-');
                }
                for(size_t s = 0; s < seqs2.size(); s++) {
                    aln_seqs[seqs1.size() + s].insert(0, 1, seqs2[s][j - 1]);
                }
                j--;
            }
            for(size_t s = 0; s < seqs1.size(); s++) {
                aln_seqs[s].insert(0, 1, '-');
            }
            for(size_t s = 0; s < seqs2.size(); s++) {
                aln_seqs[seqs1.size() + s].insert(0, 1, seqs2[s][j - 1]);
            }
            j--;
            break;
        case 2:  // deletion
            while(Bq(i, j) == 1) {
                for(size_t s = 0; s < seqs1.size(); s++) {
                    aln_seqs[s].insert(0, 1, seqs1[s][i - 1]);
                }
                for(size_t s = 0; s < seqs2.size(); s++) {
                    aln_seqs[seqs1.size() + s].insert(0, 1, '-');
                }
                i--;
            }
            for(size_t s = 0; s < seqs1.size(); s++) {
                aln_seqs[s].insert(0, 1, seqs1[s][i - 1]);
            }
            for(size_t s = 0; s < seqs2.size(); s++) {
                aln_seqs[seqs1.size() + s].insert(0, 1, '-');
            }
            i--;
            break;
        }
    }
    aln.f.seq_data.clear();
    aln.f.seq_data = aln_seqs;

    return EXIT_SUCCESS;
}

/* Weighted average of nucleotide freqs for profile */
float nuc_pi(Vector4f n, Vector5f pis) {
    float val = 0;
    for(int i = 0; i < 4; i++) {
        val += n[i] == 0 ? 0 : n[i] * exp(pis[i]);
    }

    return val;
}

TEST_CASE("[profile_aln.cc] nuc_pi") {
    Vector4f nucleotides = {0.25, 0.25, 0.25, 0.25};
    Vector5f nuc_frequencies;
    nuc_frequencies << logf(0.3), logf(0.2), logf(0.2), logf(0.3), logf(0.25);

    CHECK(nuc_pi(nucleotides, nuc_frequencies) == doctest::Approx(0.25));
    nucleotides << 1, 0, 0, 0;
    CHECK(nuc_pi(nucleotides, nuc_frequencies) == doctest::Approx(0.3));
    nucleotides << 0.5, 0.1, 0.2, 0.2;
    CHECK(nuc_pi(nucleotides, nuc_frequencies) == doctest::Approx(0.27));
}
