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

#include <coati/mutation_fst.hpp>

namespace coati {

/**
 * @brief Create Muse and Gaut codon model FST.
 *
 * @details Create an FST that represents the codon substitution model by Muse
 * and Gaut. For each codon to codon substitution value (61*61), we add three
 * arcs that model the nucleotide to nucleotide changes. The first arc carries
 * the substitution probability.
 *
 * @param[in] br_len float branch length.
 * @param[in] omega float nonsynonymous-synonymous bias.
 * @param[in] pi std::vector<coati::float_t> nucleotide frequencies (A,C,G,T).
 * @param[in] sigma std::vector<coati::float_t> transition probabilities for GTR
 * substitution model.
 *
 * @retval coati::VectorFstStdArc Muse and Gaut codon model FST.
 */
VectorFstStdArc mg94(float br_len, float omega,
                     const std::vector<coati::float_t>& pi,
                     const std::vector<coati::float_t>& sigma) {
    using coati::utils::get_nuc;
    coati::Matrixf P = mg94_p(br_len, omega, pi, sigma);

    // Add state 0 and make it the start state
    VectorFstStdArc mg94;
    mg94.AddState();
    mg94.SetStart(0);

    // Create FST
    int r = 1;
    for(uint8_t i = 0; i < 61; i++) {
        for(uint8_t j = 0; j < 61; j++) {
            add_arc(mg94, 0, r, get_nuc(i, 0) + 1, get_nuc(j, 0) + 1, P(i, j));
            add_arc(mg94, r, r + 1, get_nuc(i, 1) + 1, get_nuc(j, 1) + 1);
            add_arc(mg94, r + 1, 0, get_nuc(i, 2) + 1, get_nuc(j, 2) + 1);
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
// GCOVR_EXCL_START
TEST_CASE("mg94") {
    VectorFstStdArc mut_fst(mg94(0.0133, 0.2, {0.308, 0.185, 0.199, 0.308}));

    REQUIRE(Verify(mut_fst));          // openfst built-in sanity check
    CHECK_EQ(mut_fst.NumArcs(0), 16);  // 4x4 nuc to nuc arcs from start state
}
// GCOVR_EXCL_STOP

/**
 * @brief Create dna marginal Muse and Gaut codon model FST
 *
 * @details Given the codon substitution model by Muse and Gaut, marginalize
 * it down to a 4x4 nucleotide substitution model. Algorithm loops over all
 * codon to codon substitutions and adds the probabilities of all the times
 * a nucleotide replaces another (e.g. P(A|A) = P(AAA|AAA) + P(CCA|CCA) + ... ).
 * Probabilities are then normalized.
 *
 * @param[in] br_len float branch length.
 * @param[in] omega float nonsynonymous-synonymous bias.
 * @param[in] pi std::vector<coati::float_t> nucleotide frequencies (A,C,G,T).
 *
 * @retval coati::VectorFstStdArc dna marginal Muse and Gaut codon model FST;
 */
VectorFstStdArc dna(float br_len, float omega,
                    const std::vector<coati::float_t>& pi) {
    using coati::utils::get_nuc;

    coati::Matrixf P = mg94_p(br_len, omega, pi);

    // Add state 0 and make it the start state
    VectorFstStdArc dna;
    dna.AddState();
    dna.SetStart(0);

    coati::Matrixf dna_p(4, 4);

    for(uint8_t cod = 0; cod < 61; cod++) {  // for each codon
        for(int pos = 0; pos < 3; pos++) {   // for each position in a codon
            for(int nucf = 0; nucf < 4; nucf++) {  // for each nucleotide (from)
                for(int nuct = 0; nuct < 4; nuct++) {  // for each nuc (to)
                    for(int i = 0; i < 61; i++) {      // sum over all codons
                        dna_p(nucf, nuct) +=
                            (get_nuc(i, pos) == nuct
                                 ? get_nuc(cod, pos) == nucf ? P(cod, i) : 0.0f
                                 : 0.0f);
                    }
                }
            }
        }
    }

    float row_sums[4]{0.0f, 0.0f, 0.0f, 0.0f};
    // normalize probabilities and add an arc to the FST (one for each nuc).
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
// GCOVR_EXCL_START
TEST_CASE("dna") {
    // float branch_length = 0.0133;
    VectorFstStdArc dna_fst = dna(0.0133, 0.2, {0.308, 0.185, 0.199, 0.308});

    REQUIRE(Verify(dna_fst));          // openfst built-in sanity check
    CHECK_EQ(dna_fst.NumArcs(0), 16);  // all 4x4 nuc transitions
    CHECK_EQ(dna_fst.NumStates(), 1);  // trans are indp, only 1 state needed

    // NOLINTNEXTLINE(clang-diagnostic-unused-variable)
    float dna_val[16]{0.9961381369, 0.0005952569, 0.0028695324, 0.0003970738,
                      0.0009135811, 0.9933360211, 0.0008441978, 0.0049061999,
                      0.0042575611, 0.0008198302, 0.9941940598, 0.0007285488,
                      0.0003882735, 0.0031330203, 0.0004814705, 0.9959972357};
    fst::StateIterator<fst::StdFst> siter(dna_fst);  // FST state iterator
    fst::ArcIteratorData<fst::StdArc> data;
    dna_fst.InitArcIterator(siter.Value(), &data);

    for(auto i = 0; i < 16; i++) {
        CHECK_EQ(data.arcs[i].weight.Value(),
                 doctest::Approx(-::logf(dna_val[i])));
    }
}
// GCOVR_EXCL_STOP

/**
 * @brief Create affine gap indel model FST.
 *
 * @details Indel FST model:
 *  (note: matches 6 -> 0 omitted)
 *                  -------------------------
 *                  |     del extension     |
 *         deletion (4) <-------> (5)       |
 *                  ^              |        |
 *  start           |             \/        \/
 *   (0) --------> (3) ------->  (6) ----> (7) end
 *    |             ^
 *   \/             |
 *   (1) <-------> (2)
 * insertion    ins extension
 *
 * @param[in] gap_open float gap opening score.
 * @param[in] gap_extend float gap extension score.
 * @param[in] pi std::vector<float> nucleotide frequencies (A,C,G,T).
 *
 * @retval coati::VectorFstStdArc indel model FST.
 */
VectorFstStdArc indel(float gap_open, float gap_extend,
                      const std::vector<float>& pi, float bc_error) {
    VectorFstStdArc indel_fst;
    int start = 0, insertion = 1, deletion = 4, match = 6, end = 7;

    // Add state 0 and make it the start state
    indel_fst.AddState();
    indel_fst.SetStart(start);

    // Insertion
    add_arc(indel_fst, start, insertion, 0, 0, gap_open);  // label 0 is <eps>
    add_arc(indel_fst, start, 3, 0, 0, 1.0f - gap_open);

    for(int i = 0; i < 4; i++) {
        add_arc(indel_fst, insertion, 2, 0, i + 1, pi[i]);
    }

    add_arc(indel_fst, insertion, 2, 0, 5);  // 5 as ilabel/olabel is N
    add_arc(indel_fst, 2, insertion, 0, 0, gap_extend);
    add_arc(indel_fst, 2, 3, 0, 0, 1.0f - gap_extend);

    // Deletion
    add_arc(indel_fst, 3, deletion, 0, 0, gap_open);
    add_arc(indel_fst, 3, match, 0, 0, 1.0f - gap_open);

    for(int i = 0; i < 4; i++) {
        add_arc(indel_fst, deletion, 5, i + 1);
    }

    add_arc(indel_fst, 5, deletion, 0, 0, gap_extend);
    add_arc(indel_fst, 5, match, 0, 0, 1.0f - gap_extend);

    // Matches
    for(int i = 1; i < 5; i++) {
        add_arc(indel_fst, match, start, i, i, 1 - 3 * bc_error);  // nuc -> nuc
        add_arc(indel_fst, match, start, i, 5);                    // nuc -> N
    }

    // Base calling error
    for(int i = 1; i < 5; ++i) {
        for(int j = 1; j < 5; ++j) {
            if(i == j) {
                continue;
            }
            add_arc(indel_fst, match, start, i, j, bc_error);
        }
    }

    // End probabilities
    add_arc(indel_fst, match, end);                       // match to end
    add_arc(indel_fst, 5, end, 0, 0, 1.0f - gap_extend);  // del to end

    // Set final state & optimize
    indel_fst.SetFinal(end, 0.0);

    VectorFstStdArc indel_rmep;
    indel_rmep = fst::RmEpsilonFst<fst::StdArc>(indel_fst);  // epsilon removal

    return optimize(indel_rmep);
}

/// @private
// GCOVR_EXCL_START
TEST_CASE("indel") {
    VectorFstStdArc indel_model(
        indel(0.001, 1.f - 1.f / 6.f, {0.308, 0.185, 0.199, 0.308}, 0.0001));

    REQUIRE(Verify(indel_model));  // openfst built-in sanity check
    CHECK_EQ(indel_model.NumStates(), 4);
    CHECK_EQ(indel_model.NumArcs(0), 29);  // number of outcoming arcs
    CHECK_EQ(indel_model.NumArcs(1), 29);
    CHECK_EQ(indel_model.NumArcs(2), 24);
    CHECK_EQ(indel_model.NumArcs(3), 1);
}
// GCOVR_EXCL_STOP

/**
 * @brief Add arc to FST.
 *
 * @param[in,out] fst coati::VectorFstStdArc FST to be added an arc.
 * @param[in] src int source state.
 * @param[in] dest int destination state.
 * @param[in] ilabel int input label.
 * @param[in] olabel int output label.
 * @param[in] score float score of the arc.
 */
void add_arc(VectorFstStdArc& fst, int src, int dest, int ilabel, int olabel,
             float score) {
    if(score == 1.0) {
        score = 0.0;
    } else if(score == 0.0) {
        score = static_cast<float>(INT_MAX);
    } else {
        score = -logf(score);
    }

    if(fst.NumStates() <= dest) {
        fst.AddState();
        fst.AddArc(src, fst::StdArc(ilabel, olabel, score, dest));
    } else {
        fst.AddArc(src, fst::StdArc(ilabel, olabel, score, dest));
    }
}

/**
 * @brief Create FSA (acceptor) from a sequence.
 *
 * @param[in] content std::string sequence to be converted to an FSA.
 * @param[in,out] accept coati::VectorFstStdArc empty FSA.
 *
 * @retval true if run successfully.
 */
bool acceptor(const std::string_view content, VectorFstStdArc& accept) {
    std::map<char, int> syms = {
        {'-', 0}, {'A', 1}, {'C', 2}, {'G', 3}, {'T', 4}, {'U', 4}, {'N', 5},
        {'a', 1}, {'c', 2}, {'g', 3}, {'t', 4}, {'u', 4}, {'n', 5}};

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
 * @brief Optimize FST: remove epsilons, determinize, and minimize.
 *
 * @param[in] fst_raw coati::VectorFstStdArc FST to be optimized.
 *
 * @retval coati::VectorFstStdArc optimized FST.
 */
VectorFstStdArc optimize(VectorFstStdArc& fst_raw) {
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
