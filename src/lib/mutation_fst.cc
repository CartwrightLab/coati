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

#include <coati/mutation_fst.hpp>

using namespace std;
using namespace fst;

/* Create Muse and Gaut codon model FST */
void mg94(VectorFstStdArc& mut_fst, float br_len) {
    Matrix64f P;
    mg94_p(P, br_len);

    // Add state 0 and make it the start state
    VectorFstStdArc mg94;
    mg94.AddState();
    mg94.SetStart(0);

    // Creat FST
    int r = 1;
    for(uint8_t i = 0; i < 64; i++) {
        for(uint8_t j = 0; j < 64; j++) {
            add_arc(mg94, 0, r, ((i & 48) >> 4) + 1, ((j & 48) >> 4) + 1,
                    P(i, j));
            add_arc(mg94, r, r + 1, ((i & 12) >> 2) + 1, ((j & 12) >> 2) + 1);
            add_arc(mg94, r + 1, 0, (i & 3) + 1, (j & 3) + 1);
            r = r + 2;
        }
    }

    // Set final state & optimize
    mg94.SetFinal(0, 0.0);

    VectorFstStdArc mg94_rmep;
    mg94_rmep = RmEpsilonFst<StdArc>(mg94);  // epsilon removal

    mut_fst = optimize(mg94_rmep);
}

TEST_CASE("[mutation_fst.cc] mg94") {
    VectorFstStdArc mut_fst;
    float branch_length = 0.0133;
    mg94(mut_fst, branch_length);

    CHECK(Verify(mut_fst));           // openfst built-in sanity check
    CHECK(mut_fst.NumArcs(0) == 16);  // 4x4 nuc to nuc arcs from start state
    CHECK(mut_fst.NumStates() == 241);
}

/* Create dna marginal Muse and Gaut codon model FST*/
void dna(VectorFstStdArc& mut_fst, float br_len) {
    Matrix64f P;
    mg94_p(P, br_len);

    // Add state 0 and make it the start state
    VectorFstStdArc dna;
    dna.AddState();
    dna.SetStart(0);

    Matrix4f dna_p = Matrix4f::Zero();

    for(uint8_t cod = 0; cod < 64; cod++) {     // for each codon
        for(int pos = 0; pos < 3; pos++) {      // for each position in a codon
            for(int nuc = 0; nuc < 4; nuc++) {  // for each nucleotide (from)
                for(int nuc2 = 0; nuc2 < 4;
                    nuc2++) {                      // for each nucleotide (to)
                    for(int i = 0; i < 64; i++) {  // sum over all codons
                        dna_p(nuc, nuc2) +=
                            (((i & (uint8_t)(48 / pow(4, pos))) >>
                              (4 - 2 * pos)) == nuc2
                                 ? ((cod & (uint8_t)(48 / pow(4, pos))) >>
                                    (4 - 2 * pos)) == nuc
                                       ? P(cod, i)
                                       : 0.0f
                                 : 0.0f);
                    }
                }
            }
        }
    }

    for(int i = 0; i < 4; i++) {
        dna_p.row(i) /= dna_p.row(i).sum();
        for(int j = 0; j < 4; j++) {
            add_arc(dna, 0, 0, i + 1, j + 1, dna_p(i, j));
        }
    }

    // Set final state & optimize
    dna.SetFinal(0, 0.0);
    mut_fst = optimize(dna);
}

TEST_CASE("[mutation_fst.cc] dna") {
    VectorFstStdArc dna_fst;
    float branch_length = 0.0133;
    dna(dna_fst, branch_length);

    CHECK(Verify(dna_fst));           // openfst built-in sanity check
    CHECK(dna_fst.NumArcs(0) == 16);  // all 4x4 nuc transitions
    CHECK(dna_fst.NumStates() ==
          1);  // all transitions are indp, only 1 state needed
}

/* Create FST that maps nucleotide to AA position */
void nuc2pos(VectorFstStdArc& n2p) {
    // Add state 0 and make it the start state
    n2p.AddState();
    n2p.SetStart(0);

    int state = 1;  // variable to keep track of states
    int cod = 101;  // variable to keep track of codons

    for(int i = 1; i < 5; i++) {
        for(int j = 1; j < 5; j++) {
            for(int h = 1; h < 5; h++) {
                add_arc(n2p, 0, state, i, cod);
                add_arc(n2p, state, state + 1, j, cod + 1);
                add_arc(n2p, state + 1, 0, h, cod + 2);
                state += 2;
                cod += 3;
            }
        }
    }

    n2p.SetFinal(0, 0.0);
}

TEST_CASE("[mutation_fst.cc] nuc2pos") {
    VectorFstStdArc n2p_fst;
    nuc2pos(n2p_fst);

    CHECK(Verify(n2p_fst));             // openfst built-in sanity check
    CHECK(n2p_fst.NumArcs(0) == 64);    // one position for every aminoacid (AA)
    CHECK(n2p_fst.NumStates() == 129);  // 64 AA with 2 states each + init state
}

/* Create affine gap indel model FST*/
void indel(VectorFstStdArc& indel_model, string model) {
    float deletion = 0.001, insertion = 0.001;
    float deletion_ext = 1.0 - 1.0 / 6.0, insertion_ext = 1.0 - 1.0 / 6.0;
    float nuc_freqs[2][4] = {{0.308, 0.185, 0.199, 0.308},
                             {0.2676350, 0.2357727, 0.2539630, 0.2426323}};
    int m = model.compare("ecm") == 0 ? 1 : 0;

    VectorFstStdArc indel_fst;

    // Add state 0 and make it the start state
    indel_fst.AddState();
    indel_fst.SetStart(0);

    // Insertion
    add_arc(indel_fst, 0, 1, 0, 0, insertion);  // 0 as ilabel/olabel is <eps>
    add_arc(indel_fst, 0, 3, 0, 0, 1.0f - insertion);

    for(int i = 0; i < 4; i++) {
        add_arc(indel_fst, 1, 2, 0, i + 1, nuc_freqs[m][i]);
    }

    add_arc(indel_fst, 1, 2, 0, 5);  // 5 as ilabel/olabel is N
    add_arc(indel_fst, 2, 1, 0, 0, insertion_ext);
    add_arc(indel_fst, 2, 3, 0, 0, 1.0f - insertion_ext);

    // Deletion
    add_arc(indel_fst, 3, 4, 0, 0, deletion);
    add_arc(indel_fst, 3, 6, 0, 0, 1.0f - deletion);

    for(int i = 0; i < 4; i++) {
        add_arc(indel_fst, 4, 5, i + 1);
    }

    add_arc(indel_fst, 4, 7);

    add_arc(indel_fst, 5, 4, 0, 0, deletion_ext);
    add_arc(indel_fst, 5, 6, 0, 0, 1.0f - deletion_ext);

    // Matches
    for(int i = 0; i < 4; i++) {
        add_arc(indel_fst, 6, 0, i + 1, i + 1);
        add_arc(indel_fst, 6, 0, i + 1, 5);
    }

    add_arc(indel_fst, 6, 7);

    // Set final state & optimize
    indel_fst.SetFinal(7, 0.0);

    VectorFstStdArc indel_rmep;
    indel_rmep = RmEpsilonFst<StdArc>(indel_fst);  // epsilon removal

    indel_model = optimize(indel_rmep);
}

TEST_CASE("[mutation_fst.cc] indel") {
    VectorFstStdArc indel_model;
    string model = "m-coati";

    indel(indel_model, model);

    CHECK(Verify(indel_model));  // openfst built-in sanity check
    CHECK(indel_model.NumStates() == 4);
    CHECK(indel_model.NumArcs(0) == 17);  // number of outcoming arcs
    CHECK(indel_model.NumArcs(1) == 17);
    CHECK(indel_model.NumArcs(2) == 12);
    CHECK(indel_model.NumArcs(3) == 1);
}
