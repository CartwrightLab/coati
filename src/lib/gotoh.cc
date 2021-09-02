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
                  Matrix& P_m) {
    // P matrix for marginal Muse and Gaut codon model
    Tensor p = mg94_marginal_p(P_m);

    std::string seq_a = sequences[0];
    std::string seq_b = sequences[1];
    int m = static_cast<int>(sequences[0].length());
    int n = static_cast<int>(sequences[1].length());

    // ensure that length of first sequence (reference) is multiple of 3
    if(m % 3 != 0) {
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
    D(0, 1) = -insertion - nuc_freqs[nt4_table[pos]] - no_insertion_ext;
    P(0, 1) = -insertion - nuc_freqs[nt4_table[pos]] - no_insertion_ext;
    Bd(0, 1) = 1;
    Bp(0, 1) = 2;
    D(1, 0) = -no_insertion - deletion - no_deletion_ext;
    Q(1, 0) = -no_insertion - deletion - no_deletion_ext;
    Bd(1, 0) = 2;
    Bq(1, 0) = 2;

    // fill first row of D
    if(n + 1 >= 2) {
        for(int j = 2; j < n + 1; j++) {
            pos = seq_b[j - 1];
            D(0, j) = D(0, j - 1) - insertion_ext - nuc_freqs[nt4_table[pos]];
            P(0, j) = P(0, j - 1) - insertion_ext - nuc_freqs[nt4_table[pos]];
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

    std::string codon;
    float p1{NAN}, p2{NAN}, q1{NAN}, q2{NAN}, d{NAN};

    for(int i = 1; i < m + 1; i++) {
        codon = seq_a.substr((((i - 1) / 3) * 3), 3);  // current codon
        for(int j = 1; j < n + 1; j++) {
            // insertion
            pos = seq_b[j - 1];
            p1 = P(i, j - 1) - insertion_ext - nuc_freqs[nt4_table[pos]];
            p2 = Bd(i, j - 1) == 0
                     ? D(i, j - 1) - insertion - nuc_freqs[nt4_table[pos]] -
                           no_insertion_ext
                 : Bd(i, j - 1) == 1
                     ? D(i, j - 1) - insertion_ext - nuc_freqs[nt4_table[pos]]
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
                    logf(transition(codon, (i) % 3, seq_b[j - 1], p));
            } else if(Bd(i - 1, j - 1) == 1) {
                d = D(i - 1, j - 1) - no_deletion -
                    logf(transition(codon, (i) % 3, seq_b[j - 1], p));
            } else {
                d = D(i - 1, j - 1) -
                    logf(transition(codon, (i) % 3, seq_b[j - 1], p));
            }

            // D[i,j] = highest weight between insertion, deletion, and
            // match/mismatch. In this case, lowest (-logf(weight)) value
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

    aln.weight = D(m, n);  // weight

    // backtracking to obtain alignment
    return backtracking(Bd, Bp, Bq, seq_a, seq_b, aln);
}

/* Dynamic Programming with no frameshifts*/
int gotoh_noframeshifts(std::vector<std::string> sequences, alignment_t& aln,
                        Matrix& P_m) {
    // P matrix for marginal Muse and Gaut codon model
    Tensor p = mg94_marginal_p(P_m);

    std::string seq_a = sequences[0];
    std::string seq_b = sequences[1];
    int m = static_cast<int>(sequences[0].length());
    int n = static_cast<int>(sequences[1].length());

    // ensure that length of first sequence (reference) is multiple of 3
    if((m % 3 != 0) || (n % 3 != 0)) {
        throw std::invalid_argument(
            "The length of both sequences must be a multiple of 3.");
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

    float insertion = 0.001;
    float deletion = 0.001;
    float insertion_ext = 1.0 - (1.0 / 6.0);
    float deletion_ext = 1.0 - (1.0 / 6.0);

    float nuc_freqs[5] = {0.308, 0.185, 0.199, 0.308, 0.25};

    // DP and backtracking matrices initialization

    // fill first values on D that are independent
    D(0, 0) = 0.0;
    Bd(0, 0) = 0;
    D(0, 3) = P(0, 3) =
        -logf(insertion) -
        logf(nuc_freqs[nt4_table[static_cast<unsigned char>(seq_b[0])]]) -
        logf(nuc_freqs[nt4_table[static_cast<unsigned char>(seq_b[1])]]) -
        logf(nuc_freqs[nt4_table[static_cast<unsigned char>(seq_b[2])]]) -
        logf(1.0f - insertion_ext) - 2 * logf(insertion_ext);
    D(3, 0) = Q(3, 0) = -logf(1.0f - insertion) - logf(deletion) -
                        2 * logf(deletion_ext) - logf(1.0f - deletion_ext);
    Bd(0, 3) = 1;
    Bd(3, 0) = Bp(0, 3) = Bq(3, 0) = 2;

    // fill first row of D
    if(n + 1 >= 6) {
        for(int j = 6; j < n + 1; j += 3) {
            D(0, j) = P(0, j) =
                D(0, j - 3) - 3 * logf(insertion_ext) -
                logf(nuc_freqs[nt4_table[static_cast<unsigned char>(
                    seq_b[j - 3])]]) -
                logf(nuc_freqs[nt4_table[static_cast<unsigned char>(
                    seq_b[j - 2])]]) -
                logf(nuc_freqs[nt4_table[static_cast<unsigned char>(
                    seq_b[j - 1])]]);
            Bd(0, j) = 1;
            Bp(0, j) = 1;
        }
    }

    // fill first column of D
    if(m + 1 >= 6) {
        for(int i = 6; i < m + 1; i += 3) {
            D(i, 0) = D(i - 3, 0) - 3 * logf(deletion_ext);
            Q(i, 0) = Q(i - 3, 0) - 3 * logf(deletion_ext);
            Bd(i, 0) = 2;
            Bq(i, 0) = 1;
        }
    }

    std::string codon;
    float p1{NAN}, p2{NAN}, q1{NAN}, q2{NAN}, d{NAN};
    int temp{0};

    // Cells with only match/mismatch (1,1) & (2,2)
    codon = seq_a.substr(0, 3);
    for(int i = 1; i < 3; i++) {
        D(i, i) = D(i - 1, i - 1) - logf(1.0f - insertion) -
                  logf(1.0f - deletion) -
                  logf(transition(codon, i % 3, seq_b[i - 1], p));
        Bd(i, i) = 0;
    }

    // Second and third row/column (match/mismatch && insertion || deletion)
    for(int i = 1; i < 3; i++) {
        for(int j = i + 3; j < n + 1; j += 3) {  // rows
            codon = seq_a.substr(0, 3);
            // match/mismatch
            if(Bd(i - 1, j - 1) == 0) {
                d = D(i - 1, j - 1) - logf(1.0f - insertion) -
                    logf(1.0f - deletion) -
                    logf(transition(codon, (i) % 3, seq_b[j - 1], p));
            } else if(Bd(i - 1, j - 1) == 1) {
                d = D(i - 1, j - 1) - logf(1.0f - deletion) -
                    logf(transition(codon, (i) % 3, seq_b[j - 1], p));
            } else {
                d = D(i - 1, j - 1) -
                    logf(transition(codon, (i) % 3, seq_b[j - 1], p));
            }
            // insertion
            p1 = P(i, j - 3) - 3 * logf(insertion_ext) -
                 logf(nuc_freqs[nt4_table[static_cast<unsigned char>(
                     seq_b[j - 3])]]) -
                 logf(nuc_freqs[nt4_table[static_cast<unsigned char>(
                     seq_b[j - 2])]]) -
                 logf(nuc_freqs[nt4_table[static_cast<unsigned char>(
                     seq_b[j - 1])]]);
            p2 = Bd(i, j - 3) == 0
                     ? D(i, j - 3) - logf(insertion) - 2 * logf(insertion_ext) -
                           logf(nuc_freqs[nt4_table[static_cast<unsigned char>(
                               seq_b[j - 3])]]) -
                           logf(nuc_freqs[nt4_table[static_cast<unsigned char>(
                               seq_b[j - 2])]]) -
                           logf(nuc_freqs[nt4_table[static_cast<unsigned char>(
                               seq_b[j - 1])]]) -
                           logf(1.0f - insertion_ext)
                 : Bd(i, j - 1) == 1
                     ? D(i, j - 1) - 3 * logf(insertion_ext) -
                           logf(nuc_freqs[nt4_table[static_cast<unsigned char>(
                               seq_b[j - 3])]]) -
                           logf(nuc_freqs[nt4_table[static_cast<unsigned char>(
                               seq_b[j - 2])]]) -
                           logf(nuc_freqs[nt4_table[static_cast<unsigned char>(
                               seq_b[j - 1])]])
                     : std::numeric_limits<float>::max();
            P(i, j) = std::fmin(p1, p2);
            Bp(i, j) =
                p1 < p2
                    ? 1
                    : 2;  // 1 is insertion extension, 2 is insertion opening
            D(i, j) = P(i, j) < d ? P(i, j) : d;
            Bd(i, j) = P(i, j) < d ? 1 : 0;
        }

        for(int j = i + 3; j < m + 1; j += 3) {            // columns
            codon = seq_a.substr((((j - 1) / 3) * 3), 3);  // current codon
            // match/mismatch
            if(Bd(j - 1, i - 1) == 0) {
                d = D(j - 1, i - 1) - logf(1.0f - insertion) -
                    logf(1.0f - deletion) -
                    logf(transition(codon, (j) % 3, seq_b[i - 1], p));
            } else if(Bd(j - 1, i - 1) == 1) {
                d = D(j - 1, i - 1) - logf(1.0f - deletion) -
                    logf(transition(codon, (j) % 3, seq_b[i - 1], p));
            } else {
                d = D(j - 1, i - 1) -
                    logf(transition(codon, (j) % 3, seq_b[i - 1], p));
            }
            // deletion
            q1 = Q(j - 3, i) - 3 * logf(deletion_ext);
            q2 = Bd(j - 3, i) == 0
                     ? D(j - 3, i) - logf(1.0f - insertion) - logf(deletion) -
                           logf(1.0f - deletion_ext) - 2 * logf(deletion_ext)
                 : Bd(j - 3, i) == 1
                     ? D(j - 3, i) - logf(1.0f - deletion_ext) -
                           logf(deletion) - 2 * logf(deletion_ext)
                     : D(j - 3, i) - 3 * logf(deletion_ext);
            Q(j, i) = std::fmin(q1, q2);
            Bq(j, i) =
                q1 < q2 ? 1
                        : 2;  // 1 is deletion extension, 2 is deletion opening
            D(j, i) = Q(j, i) < d ? Q(j, i) : d;
            Bd(j, i) = Q(j, i) < d ? 1 : 2;
        }
    }

    // Cells considering all 3 events (insertion, deletion, match/mismatch)
    for(int i = 3; i < m + 1; i++) {
        codon = seq_a.substr((((i - 1) / 3) * 3), 3);  // current codon
        for(int j = 3; j < n + 1; j += 3) {
            temp = j;
            j += i % 3;
            if(j > n) continue;
            // match/mismatch
            if(Bd(i - 1, j - 1) == 0) {
                d = D(i - 1, j - 1) - logf(1.0f - insertion) -
                    logf(1.0f - deletion) -
                    logf(transition(codon, (i) % 3, seq_b[j - 1], p));
            } else if(Bd(i - 1, j - 1) == 1) {
                d = D(i - 1, j - 1) - logf(1.0f - deletion) -
                    logf(transition(codon, (i) % 3, seq_b[j - 1], p));
            } else {
                d = D(i - 1, j - 1) -
                    logf(transition(codon, (i) % 3, seq_b[j - 1], p));
            }
            // insertion
            p1 = P(i, j - 3) - 3 * logf(insertion_ext) -
                 logf(nuc_freqs[nt4_table[static_cast<unsigned char>(
                     seq_b[j - 3])]]) -
                 logf(nuc_freqs[nt4_table[static_cast<unsigned char>(
                     seq_b[j - 2])]]) -
                 logf(nuc_freqs[nt4_table[static_cast<unsigned char>(
                     seq_b[j - 1])]]);
            p2 = Bd(i, j - 3) == 0
                     ? D(i, j - 3) - logf(insertion) - 2 * logf(insertion_ext) -
                           logf(nuc_freqs[nt4_table[static_cast<unsigned char>(
                               seq_b[j - 3])]]) -
                           logf(nuc_freqs[nt4_table[static_cast<unsigned char>(
                               seq_b[j - 2])]]) -
                           logf(nuc_freqs[nt4_table[static_cast<unsigned char>(
                               seq_b[j - 1])]]) -
                           logf(1.0f - insertion_ext)
                 : Bd(i, j - 1) == 1
                     ? D(i, j - 1) - 3 * logf(insertion_ext) -
                           logf(nuc_freqs[nt4_table[static_cast<unsigned char>(
                               seq_b[j - 3])]]) -
                           logf(nuc_freqs[nt4_table[static_cast<unsigned char>(
                               seq_b[j - 2])]]) -
                           logf(nuc_freqs[nt4_table[static_cast<unsigned char>(
                               seq_b[j - 1])]])
                     : std::numeric_limits<float>::max();
            P(i, j) = std::fmin(p1, p2);
            Bp(i, j) =
                p1 < p2
                    ? 1
                    : 2;  // 1 is insertion extension, 2 is insertion opening
            // deletion
            q1 = Q(i - 3, j) - 3 * logf(deletion_ext);
            q2 = Bd(i - 3, j) == 0
                     ? D(i - 3, j) - logf(1.0f - insertion) - logf(deletion) -
                           logf(1.0f - deletion_ext) - 2 * logf(deletion_ext)
                 : Bd(i - 3, j) == 1
                     ? D(i - 3, j) - logf(1.0f - deletion_ext) -
                           logf(deletion) - 2 * logf(deletion_ext)
                     : D(i - 3, j) - 3 * logf(deletion_ext);
            Q(i, j) = std::fmin(q1, q2);
            Bq(i, j) =
                q1 < q2 ? 1
                        : 2;  // 1 is deletion extension, 2 is deletion opening
            // D[i,j] = highest weight between insertion, deletion, and
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
            j = temp;
        }
    }

    aln.weight = D(m, n);  // weight

    // backtracking to obtain alignment
    return backtracking_noframeshifts(Bd, Bp, Bq, seq_a, seq_b, aln);
}
/* Return value from marginal MG94 model p matrix for a given transition */
float transition(const std::string& codon, int position, unsigned char nuc,
                 const Tensor& p) {
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

TEST_CASE("[gotoh.cc] transition") {
    std::string codon{"AAA"};
    Matrix P(mg94_p(0.0133));

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

/* Recover alignment given backtracking matrices for DP alignment */
int backtracking_noframeshifts(const Matrix& Bd, const Matrix& Bp,
                               const Matrix& Bq, std::string seqa,
                               std::string seqb, alignment_t& aln) {
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
                for(int h = 0; h < 3; h++) {
                    aln.f.seq_data[0].insert(0, 1, '-');
                    aln.f.seq_data[1].insert(0, 1, seqb[j - 1]);
                    j--;
                }
            }
            for(int h = 0; h < 3; h++) {
                aln.f.seq_data[0].insert(0, 1, '-');
                aln.f.seq_data[1].insert(0, 1, seqb[j - 1]);
                j--;
            }
            // deletion
        } else {
            while(Bq(i, j) == 1) {
                for(int h = 0; h < 3; h++) {
                    aln.f.seq_data[0].insert(0, 1, seqa[i - 1]);
                    aln.f.seq_data[1].insert(0, 1, '-');
                    i--;
                }
            }
            for(int h = 0; h < 3; h++) {
                aln.f.seq_data[0].insert(0, 1, seqa[i - 1]);
                aln.f.seq_data[1].insert(0, 1, '-');
                i--;
            }
        }
    }

    return 0;
}
