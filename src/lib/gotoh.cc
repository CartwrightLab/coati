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

#include <coati/gotoh.hpp>

/* Dynamic Programming implementation of Marginal MG94 model*/
int mg94_marginal(std::vector<std::string> sequences, alignment_t& aln,
                  Matrix& P_m, bool frameshifts) {
    // P matrix for marginal Muse and Gaut codon model
    Tensor p = mg94_marginal_p(P_m);

    std::string seq_a = sequences[0];
    std::string seq_b = sequences[1];
    int m = static_cast<int>(sequences[0].length());
    int n = static_cast<int>(sequences[1].length());

    // ensure that length of first sequence (reference) is multiple of 3
    if((m % 3 != 0) || (!frameshifts && (n % 3 != 0))) {
        throw std::invalid_argument(
            "Reference coding sequence length must be a multiple of 3 (" +
            std::to_string(m) + "). Exiting!");
    }

    // DP matrices for match/mismatch (D), insertion (P), and deletion (Q)
    Matrix D(m + 1, n + 1, std::numeric_limits<float>::max());
    Matrix P(m + 1, n + 1, std::numeric_limits<float>::max());
    Matrix Q(m + 1, n + 1, std::numeric_limits<float>::max());

    // backtracking info matrices for match/mismatch (Bd), insert (Bp), and
    // deletion (Bq)
    Matrix Bd(m + 1, n + 1, -1.0f);
    Matrix Bp(m + 1, n + 1, -1.0f);
    Matrix Bq(m + 1, n + 1, -1.0f);

    float insertion = logf(0.001);
    float deletion = logf(0.001);
    float insertion_ext = logf(1.0f - (1.0 / 6.0));
    float deletion_ext = logf(1.0f - (1.0 / 6.0));
    float no_insertion = logf(1.0f - 0.001);
    float no_deletion = logf(1.0f - 0.001);
    float no_insertion_ext = logf(1.0f / 6.0);
    float no_deletion_ext = logf(1.0f / 6.0);

    float nuc_freqs[5] = {logf(0.308), logf(0.185), logf(0.199), logf(0.308),
                          logf(0.25)};

    // DP and backtracking matrices initialization

    // fill first values on D that are independent
    unsigned char pos = seq_b[0];
    D(0, 0) = 0.0;  // 0.0;
    Bd(0, 0) = 0;
    D(0, 1) = -insertion - nuc_freqs[nt4_table[pos]];
    Bd(0, 1) = 1;
    D(1, 0) = -no_insertion - deletion;
    Bd(1, 0) = 2;

    // fill first row of D
    if(n + 1 >= 2) {
        for(int j = 2; j < n + 1; j++) {
            pos = seq_b[j - 1];
            D(0, j) = D(0, j - 1) - insertion_ext - nuc_freqs[nt4_table[pos]];
            Bd(0, j) = 1;
        }
    }

    // fill first column of D
    if(m + 1 >= 2) {
        for(int i = 2; i < m + 1; i++) {
            D(i, 0) = D(i - 1, 0) - deletion_ext;
            Bd(i, 0) = 2;
        }
    }

    std::string codon;
    float p1{NAN}, p2{NAN}, q0{NAN}, q1{NAN}, q2{NAN}, d{NAN}, mch{NAN};

    for(int i = 1; i < m + 1; i++) {
        codon = seq_a.substr((((i - 1) / 3) * 3), 3);  // current codon
        for(int j = 1; j < n + 1; j++) {
            // insertions can follow matches or insertions
            pos = seq_b[j - 1];
            p1 = P(i, j - 1) - insertion_ext -
                 nuc_freqs[nt4_table[pos]];  // gap extend
            p2 = D(i, j - 1) - insertion -
                 nuc_freqs[nt4_table[pos]];  // gap open
            P(i, j) = std::fmin(p1, p2);
            // 1 is insertion extension, 2 is insertion opening
            Bp(i, j) = p1 < p2 ? 1 : 2;

            // deletions can follow matches, insertions, or deletions
            q0 = Q(i - 1, j) - deletion_ext;
            q1 = D(i - 1, j) - no_insertion - deletion;
            q2 = P(i - 1, j) - no_insertion_ext - deletion;
            Q(i, j) = std::fmin(std::fmin(q1, q2), q0);
            // 1 is deletion extension, 2 is deletion opening
            Bq(i, j) = ((q0 < q1) && (q0 < q2)) ? 1 : 2;

            // matches can follow matches, insertions, or deletions
            mch = logf(transition(codon, (i) % 3, seq_b[j - 1], p, j % 3,
                                  frameshifts));
            if(Bd(i - 1, j - 1) == 0) {  // match -> match
                d = D(i - 1, j - 1) - no_insertion - no_deletion - mch;
            } else if(Bd(i - 1, j - 1) == 1) {  // insertion -> match
                d = D(i - 1, j - 1) - no_insertion_ext - no_deletion - mch;
            } else {  // deletion -> match
                d = D(i - 1, j - 1) - no_deletion_ext - mch;
            }

            // D[i,j] = highest weight between insertion, deletion, and
            // match/mismatch. In this case, lowest (-logf(weight)) value
            float path = 0.f;
            auto score = d;
            if(P(i, j) < score) {
                score = P(i, j);
                path = 1;
            }
            if(Q(i, j) < score) {
                score = Q(i, j);
                path = 2;
            }
            D(i, j) = score;
            Bd(i, j) = path;
        }
    }

    // adjust terminal state
    D(m, n) -= no_insertion;
    P(m, n) -= no_insertion_ext;
    float path = 0.f;
    auto score = D(m, n);
    if(P(m, n) < score) {
        score = P(m, n);
        path = 1;
    }
    if(Q(m, n) < score) {
        score = Q(m, n);
        path = 2;
    }
    Bd(m, n) = path;
    D(m, n) = score;

    aln.weight = D(m, n);  // weight

    // backtracking to obtain alignment
    return backtracking(Bd, Bp, Bq, seq_a, seq_b, aln);
}

/* Return value from marginal MG94 model p matrix for a given transition */
float transition(const std::string& codon, int position, unsigned char nuc,
                 const Tensor& p, int position2, bool frameshifts) {
    if(!frameshifts && (position != position2)) {
        return 0;
    }

    position = position == 0 ? 2 : position == 1 ? 0 : 1;

    if(nuc != 'N') {
        return p(cod_int(codon), position, nt4_table[nuc]);
    }
    float val = 0.0;
    for(int i = 0; i < 4; i++) {
        val += p(cod_int(codon), position, i);
    }
    return val / 4.0f;
}

TEST_CASE("transition") {
    std::string codon{"AAA"};
    Matrix P(mg94_p(0.0133, 0.2));

    CHECK(transition(codon, 1, 'N', mg94_marginal_p(P)) ==
          doctest::Approx(0.25));
}

/* Recover alignment given backtracking matrices for DP alignment */
int backtracking(const Matrix& Bd, const Matrix& Bp, const Matrix& Bq,
                 std::string seqa, std::string seqb, alignment_t& aln) {
    int i = static_cast<int>(seqa.length());
    int j = static_cast<int>(seqb.length());

    // vector<string> alignment;
    aln.f.seq_data.emplace_back();
    aln.f.seq_data.emplace_back();

    while((i != 0) || (j != 0)) {
        // match/mismatch
        if(Bd(i, j) == 0) {
            aln.f.seq_data[0].insert(0, 1, seqa[i - 1]);
            aln.f.seq_data[1].insert(0, 1, seqb[j - 1]);
            i--;
            j--;
            // insertion
        } else if(Bd(i, j) == 1) {
            while(Bp(i, j) == 1) {
                aln.f.seq_data[0].insert(0, 1, '-');
                aln.f.seq_data[1].insert(0, 1, seqb[j - 1]);
                j--;
            }
            aln.f.seq_data[0].insert(0, 1, '-');
            aln.f.seq_data[1].insert(0, 1, seqb[j - 1]);
            j--;
            // deletion
        } else {
            while(Bq(i, j) == 1) {
                aln.f.seq_data[0].insert(0, 1, seqa[i - 1]);
                aln.f.seq_data[1].insert(0, 1, '-');
                i--;
            }
            aln.f.seq_data[0].insert(0, 1, seqa[i - 1]);
            aln.f.seq_data[1].insert(0, 1, '-');
            i--;
        }
    }

    return 0;
}
