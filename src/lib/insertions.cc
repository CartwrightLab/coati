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

#include <coati/insertions.hpp>

namespace coati {
/**
 * @brief Store open insertions in a sparse vector given two aligned sequences.
 *
 * @param[in] ref std::string reference/ancestor sequence.
 * @param[in] seq std::string descendant sequence.
 * @param[in, out] insertions_vector coati::SparseVectorInt vector to store open
 *  insertion positions.
 *
 * @retval true if run was successful, false otherwise.
 */
bool insertion_flags(const std::string_view ref, const std::string_view seq,
                     SparseVectorInt& insertions_vector) {
    if(ref.length() != seq.length()) {
        return false;  // return if length is diff
    }

    if(ref.find('-') == std::string::npos) {
        return true;  // return if no insertions
    }

    for(std::size_t i = 0; i < ref.length(); i++) {
        if(ref[i] == '-') {
            insertions_vector.coeffRef(static_cast<int64_t>(i)) =
                111;  // 'o' (open) in ASCII dec value
        }
    }

    return true;
}

/// @private
// GCOVR_EXCL_START
TEST_CASE("insertion_flags") {
    SparseVectorInt insertions(7);

    SUBCASE("Different length - fail") {
        REQUIRE_FALSE(insertion_flags("TCA-TC", "TCAGTCG", insertions));
    }

    SUBCASE("Two insertions") {
        REQUIRE(insertion_flags("TCA-TC-", "TCAGTCG", insertions));
        CHECK_EQ(insertions.nonZeros(), 2);
        CHECK_EQ(insertions.coeffRef(3), 'o');
        CHECK_EQ(insertions.coeffRef(6), 111);
    }
}
// GCOVR_EXCL_STOP

/**
 * @brief Collapse insertions between two sequences.
 *
 * @details Collapse insertions between two descendant sequences, respectively
 * pairwise aligned to ref sequence. If an insertion is closed, gaps are added
 * to other sequences. However, if there are multiple open insertions with same
 * nucleotide, then the insertion is closed. Unless all other sequences have
 * open insertions, in which case they remain as open.
 *
 * @param[in] ins_data coati::insertion_vector vector with sequences and
 *  insertions to be merged.
 * @param[in,out] merged_data insertion_data_t names, sequences, and positions
 *  of merged insertions.
 */
void merge_indels(coati::insertion_vector& ins_data,
                  insertion_data_t& merged_data) {
    // SparseVectorInt is expected to be of 2*length(seq) so that
    //   after insertions being merged there is no index out of range/bounds

    if(ins_data.size() < 2) {
        throw std::runtime_error("Merging indels of only 1 sequence.");
    }

    // count total number of insertion flags to process
    uint64_t num_gaps = 0, processed_gaps = 0;
    for(auto in_d : ins_data) {
        in_d.insertions.prune(1);
        num_gaps += static_cast<uint64_t>(in_d.insertions.nonZeros());
    }

    std::size_t pos = 0;
    while(processed_gaps < num_gaps) {  // while gaps not processed
        // (i) for every closed insertion --> add gap to other sequences
        processed_gaps += add_closed_ins(ins_data, pos);

        // (ii) if insertions at every seq & all open & same char: continue
        if(check_all_open(ins_data, pos)) {
            pos++;
            processed_gaps += ins_data.size();
            continue;
        }

        // (iii) for every open ins --> find all open ins w/ same char,
        std::vector<std::size_t> indexes = find_open_ins(ins_data, pos);
        if(indexes.size() > 0) {
            add_gap(ins_data, indexes, pos);
            processed_gaps += indexes.size();
        }

        pos++;
    }
    // populate output variable with results
    // first sequences and names
    for(auto dat : ins_data) {
        for(std::size_t i = 0; i < dat.sequences.size(); i++) {
            merged_data.sequences.push_back(dat.sequences[i]);
            merged_data.names.push_back(dat.names[i]);
        }
    }
    // insertions vector should be the same for all of them, so pick one
    merged_data.insertions = ins_data[0].insertions;
}

/**
 * @brief Add closed insertions to other sequences.
 *
 * @param[in,out] ins_data coati::insertion_vector insertion data.
 * @param[in] pos std::size_t current position in alignment/sequence.
 *
 * @retval uint64_t number of gaps processed.
 */
uint64_t add_closed_ins(coati::insertion_vector& ins_data, std::size_t pos) {
    uint64_t processed_gaps{0};
    // foreach set of seqs
    for(std::size_t seq = 0; seq < ins_data.size(); seq++) {
        // if insertion is closed add a gap to the rest of sequences
        if(ins_data[seq].insertions.coeffRef(static_cast<int64_t>(pos)) == 99) {
            add_gap(ins_data, {seq}, pos);  // add gaps to rest of seqs
            seq--;  // start over checking for closed gaps
            pos++;  // move to next pos to not count ins added in this iteration
            processed_gaps++;
        }
    }
    return processed_gaps;
}

/**
 * @brief Check if all sequences have an open insertion at a given position.
 *
 * @details Check if all sequences at a given position have an open insertion
 * and the inserted nucleotide is the same.
 *
 * @param[in,out] ins_data coati::insertion_vector insertion data.
 * @param[in] pos std::size_t current position in alignment/sequence.
 *
 * @retval true if run is successful, false otherwise.
 */
bool check_all_open(coati::insertion_vector& ins_data, std::size_t pos) {
    char nuc = '0';
    for(auto& seq : ins_data) {  // foreach set of sequences
        if(pos > seq.sequences[0].length()) {
            // if position is out of bounds return false
            return false;
        }
        if(nuc == '0') {
            // on the first check assign nucleotide
            nuc = seq.sequences[0][pos];
        }
        if((seq.insertions.coeffRef(static_cast<int64_t>(pos)) != 111) ||
           (seq.sequences[0][pos] != nuc)) {
            // if position is not open or nucleotide is different, return false
            return false;
        }
    }
    return true;
}

/**
 * @brief Find open insertions with same nucleotide at a given position.
 *
 * @param[in,out] ins_data coati::insertion_vector insertion data.
 * @param[in] pos std::size_t current position in alignment/sequence.
 *
 * @retval std::vector<std::size_t> indexes of sequences with open insertion
 * and same nucleotide.
 */
std::vector<std::size_t> find_open_ins(coati::insertion_vector& ins_data,
                                       std::size_t pos) {
    std::vector<std::size_t> indexes;
    char nuc = '0';
    // foreach set of seqs
    for(std::size_t seq = 0; seq < ins_data.size(); seq++) {
        if(ins_data[seq].insertions.coeffRef(static_cast<int64_t>(pos)) ==
           111) {  // if position is open
            if(pos > ins_data[seq].sequences[0].length()) {
                // if position is over length of sequence skip
                continue;
            }
            if(nuc == '0') {
                // first time assign nucleotide and add sequence index
                nuc = ins_data[seq].sequences[0][pos];
                indexes.push_back(seq);
            } else {
                if(ins_data[seq].sequences[0][pos] == nuc) {
                    // if nucleotide is the same, add to list
                    indexes.push_back(seq);
                }
            }
        }
    }
    return indexes;
}

/// @private
// GCOVR_EXCL_START
TEST_CASE("merge_indels") {
    SUBCASE("Two sequences, two insertion vectors, two merges") {
        SparseVectorInt insA(14);  // size (upperbound) should be 2*length
        SparseVectorInt insB(14);
        insA.insert(5) = 111;  // open gap at position 5
        insB.insert(3) = 111;  // open gap at position 3
        insB.insert(6) = 111;  // open gap at position 6

        insertion_data_t results;
        insertion_data_t dataA("TCATCG", "A", insA);
        insertion_data_t dataB("TCAGTCG", "B", insB);

        coati::insertion_vector ins_data{dataA, dataB};

        merge_indels(ins_data, results);

        CHECK_EQ(results.sequences[0], "TCA-TCG");
        CHECK_EQ(results.sequences[1], "TCAGTCG");
        CHECK_EQ(results.insertions.coeffRef(3), 99);   // close gap at 3
        CHECK_EQ(results.insertions.coeffRef(6), 111);  // open gap at 6
    }

    SUBCASE("Four sequences, two insertion vectors, two merges") {
        SparseVectorInt insABC(14);  // size (upperbound) should be 2*length
        SparseVectorInt insD(14);

        insABC.insert(3) = 99;   // closed insertion in position 3
        insABC.insert(6) = 111;  // open insertion in position 6
        insD.insert(3) = 111;    // open insertion in position 3
        insD.insert(6) = 111;    // open insertion in position 6

        insertion_data_t results;
        insertion_data_t dataABC({"TCA-TCG", "TCAGTCG", "T-A-TCG"},
                                 {"A", "B", "C"}, insABC);
        insertion_data_t dataD("TCACTCG", "D", insD);

        coati::insertion_vector ins_data{dataABC, dataD};

        merge_indels(ins_data, results);

        CHECK_EQ(results.sequences[0], "TCA--TCG");
        CHECK_EQ(results.sequences[1], "TCAG-TCG");
        CHECK_EQ(results.sequences[2], "T-A--TCG");
        CHECK_EQ(results.sequences[3], "TCA-CTCG");
        CHECK_EQ(results.insertions.coeffRef(3), 99);
        CHECK_EQ(results.insertions.coeffRef(4), 99);
        CHECK_EQ(results.insertions.coeffRef(7), 111);
    }

    SUBCASE("Four sequences, two insertion vectors, three merges") {
        SparseVectorInt insABC(14);  // size (upperbound) should be 2*length
        SparseVectorInt insD(14);

        insABC.insert(3) = 99;   // closed insertion in position 3
        insABC.insert(6) = 111;  // open insertion in position 6
        insD.insert(3) = 111;    // open insertion in position 3
        insD.insert(6) = 111;    // open insertion in position 6

        insertion_data_t results;
        insertion_data_t dataABC({"TCA-TCG", "TCAGTCG", "T-A-TCG"},
                                 {"A", "B", "C"}, insABC);
        insertion_data_t dataD("TCACTCG", "D", insD);

        coati::insertion_vector ins_data{dataABC, dataD};

        merge_indels(ins_data, results);

        CHECK_EQ(results.sequences[0], "TCA--TCG");
        CHECK_EQ(results.sequences[1], "TCAG-TCG");
        CHECK_EQ(results.sequences[2], "T-A--TCG");
        CHECK_EQ(results.sequences[3], "TCA-CTCG");
        CHECK_EQ(results.insertions.coeffRef(3), 99);
        CHECK_EQ(results.insertions.coeffRef(4), 99);
        CHECK_EQ(results.insertions.coeffRef(7), 111);
    }

    SUBCASE("Two sequences, two insertion vectors, one merge") {
        SparseVectorInt insA(14);  // size (upperbound) should be 2*length
        SparseVectorInt insB(14);
        insA.insert(3) = 111;  // open gap at position 3
        insB.insert(3) = 111;  // open gap at position 3

        insertion_data_t results;
        insertion_data_t dataA("TCACTCG", "A", insA);
        insertion_data_t dataB("TCAGTCG", "B", insB);

        coati::insertion_vector ins_data{dataA, dataB};

        merge_indels(ins_data, results);

        CHECK_EQ(results.sequences[0], "TCAC-TCG");
        CHECK_EQ(results.sequences[1], "TCA-GTCG");
        CHECK_EQ(results.insertions.coeffRef(3), 99);  // close gap at 3
        CHECK_EQ(results.insertions.coeffRef(4), 99);  // open gap at 4
    }

    SUBCASE("Three sequences, one insertion (len 21) in C sequence") {
        SparseVectorInt insH(54);  // human
        SparseVectorInt insG(54);  // gorilla
        SparseVectorInt insC(96);  // chimp

        for(int i = 27; i < 48; i++) {
            insC.insert(i) = 111;  // open gap
        }

        insertion_data_t results;
        insertion_data_t dataH("AAATTCCAACAACATAAACAAATCTGA", "H", insH);
        insertion_data_t dataG("AAATTCCAACAACATAAACAAATCTGA", "G", insG);
        insertion_data_t dataC(
            "AAATTCCAACAACATAAACAGATCGGAAGAGAAACTATGCTTTTCTAG", "C", insC);

        coati::insertion_vector ins_data{dataH, dataG, dataC};

        merge_indels(ins_data, results);

        CHECK_EQ(results.sequences[0],
                 "AAATTCCAACAACATAAACAAATCTGA---------------------");
        CHECK_EQ(results.sequences[1],
                 "AAATTCCAACAACATAAACAAATCTGA---------------------");
        CHECK_EQ(results.sequences[2],
                 "AAATTCCAACAACATAAACAGATCGGAAGAGAAACTATGCTTTTCTAG");

        for(int i = 27; i < 48; i++) {
            CHECK_EQ(results.insertions.coeffRef(i), 99);  // closed gap
        }
    }

    SUBCASE("Two sequences, one insertion (len 3), one open insertion") {
        SparseVectorInt insA(34);
        SparseVectorInt insB(34);

        insA.insert(4) = 111;  // open gap
        insB.insert(2) = 111;
        insB.insert(3) = 111;
        insB.insert(4) = 111;
        insB.insert(7) = 111;  // same insertion as insA[4]

        insertion_data_t results;
        insertion_data_t dataA("CTTGCAT", "A", insA);
        insertion_data_t dataB("CTACGTGCAT", "B", insB);
        coati::insertion_vector ins_data{dataA, dataB};

        merge_indels(ins_data, results);

        CHECK_EQ(results.sequences[0], "CT---TGCAT");
        CHECK_EQ(results.sequences[1], "CTACGTGCAT");
        CHECK_EQ(results.insertions.coeffRef(2), 99);   // closed gap
        CHECK_EQ(results.insertions.coeffRef(3), 99);   // closed gap
        CHECK_EQ(results.insertions.coeffRef(4), 99);   // closed gap
        CHECK_EQ(results.insertions.coeffRef(7), 111);  // open gap
    }

    SUBCASE("One sequences - fail") {
        SparseVectorInt insA(34);

        insA.insert(4) = 111;  // open gap

        insertion_data_t results;
        insertion_data_t dataA("CTTGCAT", "A", insA);
        coati::insertion_vector ins_data{dataA};

        CHECK_THROWS_AS(merge_indels(ins_data, results), std::runtime_error);
    }
}
// GCOVR_EXCL_STOP

/**
 * @brief Add closed gaps to sequences and insertion vectors.
 *
 * @param[in,out] ins_data coati::insertion_vector sequences and insertion
 *  positions and type.
 * @param[in] seq_indexes std::vector<int> indexes of which sequences to add
 *  closed insertion.
 * @param[in] pos std::size_t position in sequence where to add closed
 * insertion.
 */
void add_gap(coati::insertion_vector& ins_data,
             const std::vector<std::size_t>& seq_indexes, std::size_t pos) {
    std::vector<int> add;
    std::vector<int> aux(ins_data.size());
    iota(begin(aux), end(aux), 0);

    std::set_difference(aux.begin(), aux.end(), seq_indexes.begin(),
                        seq_indexes.end(), inserter(add, add.begin()));

    // make sure insertion is closed
    for(auto seq : seq_indexes) {
        ins_data[seq].insertions.coeffRef(static_cast<int64_t>(pos)) = 99;
    }

    // add closed insertion
    for(auto seq : add) {
        for(auto& sequence : ins_data[seq].sequences) {
            // add gap to nucleotide sequence
            sequence.insert(sequence.begin() + static_cast<int64_t>(pos), '-');
        }
        // add closed gap to insertion SparseVector
        for(Eigen::Index i = ins_data[seq].insertions.cols() - 1;
            i > static_cast<int64_t>(pos); i--) {
            ins_data[seq].insertions.coeffRef(i) =
                ins_data[seq].insertions.coeffRef(i - 1);
        }
        ins_data[seq].insertions.coeffRef(static_cast<int64_t>(pos)) = 99;
    }
}

/// @private
// GCOVR_EXCL_START
TEST_CASE("add_gap") {
    SUBCASE("Three sequences with three insertion vectors") {
        SparseVectorInt insA(14);  // size (upperbound) should be 2*length
        SparseVectorInt insB(14);
        SparseVectorInt insC(14);
        insC.insert(1) = 111;

        insertion_data_t dataA("TCATCG", "A", insA);
        insertion_data_t dataB("TCAGTCG", "B", insB);
        insertion_data_t dataC("TTCATCG", "C", insC);

        std::vector<std::size_t> indexes = {2};

        coati::insertion_vector ins_data = {dataA, dataB, dataC};

        add_gap(ins_data, indexes, 1);
        CHECK_EQ(ins_data[0].sequences[0], "T-CATCG");
        CHECK_EQ(ins_data[0].insertions.coeffRef(1), 99);
        CHECK_EQ(ins_data[1].sequences[0], "T-CAGTCG");
        CHECK_EQ(ins_data[1].insertions.coeffRef(1), 99);
        CHECK_EQ(ins_data[2].sequences[0], "TTCATCG");
        CHECK_EQ(ins_data[2].insertions.coeffRef(1), 99);
    }

    SUBCASE("Two sequences with two insertion vectors - end boundaries") {
        SparseVectorInt insA(14);  // size (upperbound) should be 2*length
        SparseVectorInt insB(14);
        insB.insert(6) = 111;  // position of gap

        insertion_data_t dataA("TCATCG", "A", insA);
        insertion_data_t dataB("TCATCGC", "B", insB);

        // sequence that already have insertion
        std::vector<std::size_t> indexes = {1};

        coati::insertion_vector ins_data = {dataA, dataB};

        add_gap(ins_data, indexes, 6);

        CHECK_EQ(ins_data[0].sequences[0], "TCATCG-");
        CHECK_EQ(ins_data[0].insertions.coeffRef(6), 99);
        CHECK_EQ(ins_data[1].sequences[0], "TCATCGC");
        CHECK_EQ(ins_data[1].insertions.coeffRef(6), 99);
    }

    SUBCASE("Four sequences with 2 insertion vectors") {
        SparseVectorInt insABC(14);  // size (upperbound) should be 2*length
        SparseVectorInt insD(14);

        insABC.insert(3) = 99;   // closed insertion in position 3
        insABC.insert(6) = 111;  // open insertion in position 6
        insD.insert(3) = 111;    // open insertion in position 3
        insD.insert(6) = 111;    // open insertion in position 6

        insertion_data_t dataABC({"TCA-TCG", "TCAGTCG", "T-A-TCG"},
                                 {"A", "B", "C"}, insABC);
        insertion_data_t dataD("TCACTCG", "D", insD);

        std::vector<std::size_t> indexes = {0};

        coati::insertion_vector ins_data = {dataABC, dataD};

        add_gap(ins_data, indexes, 3);
        CHECK_EQ(ins_data[0].sequences[0], "TCA-TCG");
        CHECK_EQ(ins_data[0].sequences[1], "TCAGTCG");
        CHECK_EQ(ins_data[0].sequences[2], "T-A-TCG");
        CHECK_EQ(ins_data[0].insertions.coeffRef(3), 99);
        CHECK_EQ(ins_data[1].sequences[0], "TCA-CTCG");
        CHECK_EQ(ins_data[1].insertions.coeffRef(3), 99);
        CHECK_EQ(ins_data[1].insertions.coeffRef(4), 111);
        CHECK_EQ(ins_data[1].insertions.coeffRef(7), 111);
    }

    SUBCASE("Two sequences, one close insertion (len 3), one open insertion") {
        SparseVectorInt insA(20);
        SparseVectorInt insB(20);

        insA.insert(2) = 99;
        insA.insert(3) = 99;
        insA.insert(6) = 111;  // open gap
        insB.insert(2) = 99;
        insB.insert(3) = 99;
        insB.insert(4) = 111;
        insB.insert(7) = 111;  // same insertion as insA[4]

        insertion_data_t dataA("CT--TGCAT", "A", insA);
        insertion_data_t dataB("CTACGTGCAT", "B", insB);
        coati::insertion_vector ins_data = {dataA, dataB};
        std::vector<std::size_t> indexes = {1};

        add_gap(ins_data, indexes, 4);

        CHECK_EQ(ins_data[0].sequences[0], "CT---TGCAT");
        CHECK_EQ(ins_data[1].sequences[0], "CTACGTGCAT");
        CHECK_EQ(ins_data[0].insertions.coeffRef(0), 0);
        CHECK_EQ(ins_data[0].insertions.coeffRef(1), 0);
        CHECK_EQ(ins_data[0].insertions.coeffRef(2), 99);  // closed gap
        CHECK_EQ(ins_data[0].insertions.coeffRef(3), 99);
        CHECK_EQ(ins_data[0].insertions.coeffRef(4), 99);
        CHECK_EQ(ins_data[0].insertions.coeffRef(5), 0);  // open gap
        CHECK_EQ(ins_data[0].insertions.coeffRef(6), 0);
        CHECK_EQ(ins_data[0].insertions.coeffRef(7), 111);
        CHECK_EQ(ins_data[1].insertions.coeffRef(0), 0);
        CHECK_EQ(ins_data[1].insertions.coeffRef(1), 0);
        CHECK_EQ(ins_data[1].insertions.coeffRef(2), 99);  // closed gap
        CHECK_EQ(ins_data[1].insertions.coeffRef(3), 99);  // open gap
        CHECK_EQ(ins_data[1].insertions.coeffRef(4), 99);  // open gap
        CHECK_EQ(ins_data[1].insertions.coeffRef(5), 0);
        CHECK_EQ(ins_data[1].insertions.coeffRef(6), 0);
        CHECK_EQ(ins_data[1].insertions.coeffRef(7), 111);  // open gap
        CHECK_EQ(ins_data[1].insertions.coeffRef(8), 0);
        CHECK_EQ(ins_data[1].insertions.coeffRef(9), 0);
    }
}
// GCOVR_EXCL_STOP
}  // namespace coati
