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

namespace coati {

/**
 * \brief Create Muse and Gaut codon model FST
 *
 * @param[in] br_len float branch length.
 * @param[in] omega float nonsynonymous-synonymous bias.
 * @param[in] pi std::vector<coati::float_t> nucleotide frequencies (A,C,G,T).
 *
 * \return Muse and Gaut codon model FST (coati::VectorFstStdArc).
 */
VectorFstStdArc mg94(float br_len, float omega,
                     const std::vector<coati::float_t>& pi,
                     const std::vector<coati::float_t>& sigma) {
    coati::Matrixf P = mg94_p(br_len, omega, pi, sigma);

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
    mg94_rmep = fst::RmEpsilonFst<fst::StdArc>(mg94);  // epsilon removal

    return optimize(mg94_rmep);
}

/// @private
TEST_CASE("mg94") {
    // float branch_length = 0.0133;
    VectorFstStdArc mut_fst(mg94(0.0133, 0.2, {0.308, 0.185, 0.199, 0.308}));

    CHECK(Verify(mut_fst));           // openfst built-in sanity check
    CHECK(mut_fst.NumArcs(0) == 16);  // 4x4 nuc to nuc arcs from start state
    CHECK(mut_fst.NumStates() == 241);
}

/**
 * \brief Create dna marginal Muse and Gaut codon model FST
 *
 * @param[in] br_len float branch length.
 * @param[in] omega float nonsynonymous-synonymous bias.
 * @param[in] pi std::vector<coati::float_t> nucleotide frequencies (A,C,G,T).
 *
 * \return dna marginal Muse and Gaut codon model FST (coati::VectorFstStdArc).
 */
VectorFstStdArc dna(float br_len, float omega,
                    const std::vector<coati::float_t>& pi) {
    coati::Matrixf P = mg94_p(br_len, omega, pi);

    // Add state 0 and make it the start state
    VectorFstStdArc dna;
    dna.AddState();
    dna.SetStart(0);

    coati::Matrixf dna_p(4, 4);

    for(uint8_t cod = 0; cod < 64; cod++) {     // for each codon
        for(int pos = 0; pos < 3; pos++) {      // for each position in a codon
            for(int nuc = 0; nuc < 4; nuc++) {  // for each nucleotide (from)
                for(int nuc2 = 0; nuc2 < 4;
                    nuc2++) {                      // for each nucleotide (to)
                    for(int i = 0; i < 64; i++) {  // sum over all codons
                        dna_p(nuc, nuc2) +=
                            (((i & static_cast<uint8_t>(48 / pow(4, pos))) >>
                              (4 - 2 * pos)) == nuc2
                                 ? ((cod &
                                     static_cast<uint8_t>(48 / pow(4, pos))) >>
                                    (4 - 2 * pos)) == nuc
                                       ? P(cod, i)
                                       : 0.0f
                                 : 0.0f);
                    }
                }
            }
        }
    }

    float row_sums[4]{0.0f, 0.0f, 0.0f, 0.0f};
    for(auto i = 0; i < 4; i++) {
        for(auto j = 0; j < 4; j++) {
            row_sums[i] += dna_p(i, j);
        }
        for(auto j = 0; j < 4; j++) {
            dna_p(i, j) /= row_sums[i];
            add_arc(dna, 0, 0, i + 1, j + 1, dna_p(i, j));
        }
    }

    // Set final state & optimize
    dna.SetFinal(0, 0.0);
    return optimize(dna);
}

/// @private
TEST_CASE("dna") {
    // float branch_length = 0.0133;
    VectorFstStdArc dna_fst = dna(0.0133, 0.2, {0.308, 0.185, 0.199, 0.308});

    CHECK(Verify(dna_fst));           // openfst built-in sanity check
    CHECK(dna_fst.NumArcs(0) == 16);  // all 4x4 nuc transitions
    CHECK(dna_fst.NumStates() ==
          1);  // all transitions are indp, only 1 state needed

    // NOLINTNEXTLINE(clang-diagnostic-unused-variable)
    float dna_val[16]{
        0.0035908299, 7.56063986,    5.91386986,    7.92348003,
        7.04520988,   0.00627565011, 7.14309978,    5.38297987,
        5.47493982,   7.21499014,    0.00562130986, 7.29403019,
        7.92336988,   5.89603996,    7.7301898,     0.00355818006};
    fst::StateIterator<fst::StdFst> siter(dna_fst);  // FST state iterator
    fst::ArcIteratorData<fst::StdArc> data;
    dna_fst.InitArcIterator(siter.Value(), &data);

    for(auto i = 0; i < 16; i++) {
        CHECK(data.arcs[i].weight.Value() == doctest::Approx(dna_val[i]));
    }
}

/**
 * \brief Create FST that maps nucleotide to AA position
 *
 * \return FST to map nucleotides to AA position (coati::VectorFstStdArc).
 */
VectorFstStdArc nuc2pos() {
    // Add state 0 and make it the start state
    VectorFstStdArc n2p;
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
    return n2p;
}

/// @private
TEST_CASE("nuc2pos") {
    VectorFstStdArc n2p_fst(nuc2pos());

    CHECK(fst::Verify(n2p_fst));        // openfst built-in sanity check
    CHECK(n2p_fst.NumArcs(0) == 64);    // one position for every aminoacid (AA)
    CHECK(n2p_fst.NumStates() == 129);  // 64 AA with 2 states each + init state
}

/**
 * \brief Create affine gap indel model FST.
 *
 * @param[in] gap_open float gap opening score.
 * @param[in] gap_extend float gap extension score.
 * @param[in] pi std::vector<float> nucleotide frequencies (A,C,G,T).
 *
 * \return indel model FST (coati::VectorFstStdArc).
 */
VectorFstStdArc indel(float gap_open, float gap_extend,
                      const std::vector<float>& pi) {
    float deletion = gap_open, insertion = gap_open;
    float deletion_ext = gap_extend, insertion_ext = gap_extend;
    float nuc_freqs[4] = {pi[0], pi[1], pi[2], pi[3]};

    VectorFstStdArc indel_fst;

    // Add state 0 and make it the start state
    indel_fst.AddState();
    indel_fst.SetStart(0);

    // Insertion
    add_arc(indel_fst, 0, 1, 0, 0, insertion);  // 0 as ilabel/olabel is <eps>
    add_arc(indel_fst, 0, 3, 0, 0, 1.0f - insertion);

    for(int i = 0; i < 4; i++) {
        add_arc(indel_fst, 1, 2, 0, i + 1, nuc_freqs[i]);
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
    indel_rmep = fst::RmEpsilonFst<fst::StdArc>(indel_fst);  // epsilon removal

    return optimize(indel_rmep);
}

/// @private
TEST_CASE("indel") {
    std::string model = "m-coati";
    VectorFstStdArc indel_model(
        indel(0.001, 1.f - 1.f / 6.f, {0.308, 0.185, 0.199, 0.308}));

    CHECK(Verify(indel_model));  // openfst built-in sanity check
    CHECK(indel_model.NumStates() == 4);
    CHECK(indel_model.NumArcs(0) == 17);  // number of outcoming arcs
    CHECK(indel_model.NumArcs(1) == 17);
    CHECK(indel_model.NumArcs(2) == 12);
    CHECK(indel_model.NumArcs(3) == 1);
}

/**
 * \brief Add arc to FST
 *
 * @param[in,out] fst coati::VectorFstStdArc FST to be added an arc.
 * @param[in] src int source state.
 * @param[in] dest int destination state.
 * @param[in] ilabel int input label.
 * @param[in] olabel int output label.
 * @param[in] weight float weight of the arc.
 */
void add_arc(VectorFstStdArc& fst, int src, int dest, int ilabel, int olabel,
             float weight) {
    if(weight == 1.0) {
        weight = 0.0;
    } else if(weight == 0.0) {
        weight = static_cast<float>(INT_MAX);
    } else {
        weight = -logf(weight);
    }

    if(fst.NumStates() <= dest) {
        fst.AddState();
        fst.AddArc(src, fst::StdArc(ilabel, olabel, weight, dest));
    } else {
        fst.AddArc(src, fst::StdArc(ilabel, olabel, weight, dest));
    }
}

/**
 * \brief Create FSA (acceptor) from a fasta file
 *
 * @param[in] content std::string sequence to be converted to an FSA.
 * @param[in,out] accept coati::VectorFstStdArc empty FSA.
 */
bool acceptor(std::string content, VectorFstStdArc& accept) {
    std::map<char, int> syms = {{'-', 0}, {'A', 1}, {'C', 2}, {'G', 3},
                                {'T', 4}, {'U', 4}, {'N', 5}};

    // Add initial state
    accept.AddState();
    accept.SetStart(0);

    for(int i = 0, content_len = static_cast<int>(content.length());
        i < content_len; i++) {
        add_arc(accept, i, i + 1, syms.at(content[i]), syms.at(content[i]));
    }

    // Add final state and run an FST sanity check (Verify)
    accept.SetFinal(accept.NumStates() - 1, 0.0);
    return Verify(accept);
}

/**
 * \brief Optimize FST: remove epsilons, determinize, and minimize
 *
 * @param[in] fst_raw coati::VectorFstStdArc FST to be optimized.
 * \return optimized FST (coati::VectorFstStdArc).
 */
VectorFstStdArc optimize(VectorFstStdArc fst_raw) {
    using fst::StdArc;
    // encode FST
    fst::SymbolTable syms;
    fill_symbol_table(syms);

    fst::EncodeMapper<StdArc> encoder(fst::kEncodeLabels, fst::ENCODE);
    encoder.SetInputSymbols(&syms);
    encoder.SetOutputSymbols(&syms);
    fst::Encode(&fst_raw, &encoder);

    // reduce to more efficient form
    // 1. epsilon removal
    fst::RmEpsilon(&fst_raw);

    // 2. determinize
    VectorFstStdArc fst_det;
    fst::Determinize(fst_raw, &fst_det);

    // 3. minimize
    fst::Minimize(&fst_det);

    // decode
    fst::Decode(&fst_det, encoder);

    return fst_det;
}
}  // namespace coati
