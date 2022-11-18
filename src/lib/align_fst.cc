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

#include <coati/align_fst.hpp>
#include <filesystem>

namespace coati {

/**
 * @brief Pairwise alignment using FST composition.
 *
 * Convert input sequences as finite-state acceptors (FSA).
 * Compose seq1 FST with evolution FST then with seq2 FST. Store in aln FST.
 * Get alignment by finding shortest path on aln FST.
 *
 * @param[in] aln coati::alignment_t alignment information.
 *
 * @retval true OK
 */
bool fst_alignment(coati::alignment_t& aln) {
    using fst::StdArc;

    // scoring only works with marginal models
    if(aln.score) {
        throw std::invalid_argument("Scoring only works with marginal models.");
    }

    // read input data
    aln.data = coati::io::read_input(aln);
    if(aln.data.size() != 2) {
        throw std::invalid_argument("Exactly two sequences required.");
    }

    if(aln.seq(0).length() % 3 != 0) {
        throw std::invalid_argument(
            "Length of reference sequence must be multiple of 3.");
    }

    // set substitution matrix according to model
    coati::utils::set_subst(aln);

    // create coati FST - combines mutation and indel models
    VectorFstStdArc evolution = evo_fst(aln);

    // find alignment graph
    // 1. compose in_tape and coati FSTs
    fst::ComposeFst<StdArc> aln_inter =
        fst::ComposeFst<StdArc>(aln.data.fsts[0], evolution);
    // 2. sort intermediate composition
    VectorFstStdArc aln_inter_sort;
    aln_inter_sort = fst::ArcSortFst<StdArc, fst::OLabelCompare<StdArc>>(
        aln_inter, fst::OLabelCompare<StdArc>());
    // 3. compose intermediate and out_tape FSTs
    VectorFstStdArc graph_fst;
    fst::Compose(aln_inter_sort, aln.data.fsts[1], &graph_fst);

    // Options to find shortest path - option 3 is best:
    //
    //  case 1: ComposeFst is time: O(v1 v2 d1 (log d2 + m2)), space O(v1 v2)
    //          ShortestPath with ComposeFst (PDT) is time: O((V+E)^4),
    //          space: O((V+E)^3)
    //          Then convert to VectorFst (FST) is time, space: O(e^(O(V+E)))
    //  case 2: ComposeFst is time: O(v1 v2 d1 (log d2 + m2)), space: O(v1 v2)
    //          Convert ComposeFst (PDT) to VectorFst(FST) is time,
    //          space: O(e^(O(V+E)))
    //          Then ShortestPath with VectorFst is O(V log(V) + E)
    //  case 3: Compose is time: O(V1 V2 D1 (log D2 + M2)), space O(V1 V2 D1 M2)
    //          Then ShortestPath with VectorFst is O(V log(V) + E)

    // find shortest path through graph
    VectorFstStdArc aln_path;
    fst::ShortestPath(graph_fst, &aln_path);

    // shortestdistance = weight of shortestpath
    if(!aln.weight_file.empty()) {
        std::vector<StdArc::Weight> distance;
        std::ofstream out_w;

        fst::ShortestDistance(aln_path, &distance);
        // append weight and fasta file name info in file
        out_w.open(aln.weight_file, std::ios::app | std::ios::out);
        out_w << aln.data.path.string() << "," << aln.model << ","
              << distance[0] << std::endl;
        out_w.close();
    }

    // topsort path FST
    fst::TopSort(&aln_path);

    // write alignment
    coati::io::write_output(aln.data, aln_path);
    return true;
}

/**
 * @brief Create evolution FST - combines mutation and indel models.
 *
 * Create substitution and indel FST models.
 * Compose and optimize them to create coati FST.
 *
 * @param[in] aln coati::alignment_t raw substitution FST, nuc frequencies, and
 * gap open/extend scores.
 *
 * @return evolution FST
 */
VectorFstStdArc evo_fst(const coati::alignment_t& aln) {
    using fst::StdArc;
    // get indel FST
    VectorFstStdArc indel_fst = indel(aln.gap.open, aln.gap.extend, aln.pi);

    // sort mutation and indel FSTs
    VectorFstStdArc mutation_sort, indel_sort;
    mutation_sort = fst::ArcSortFst<StdArc, fst::OLabelCompare<StdArc>>(
        aln.subst_fst, fst::OLabelCompare<StdArc>());
    indel_sort = fst::ArcSortFst<StdArc, fst::ILabelCompare<StdArc>>(
        indel_fst, fst::ILabelCompare<StdArc>());

    // compose mutation and indel FSTs
    fst::ComposeFst<StdArc> evo_comp =
        fst::ComposeFst<StdArc>(mutation_sort, indel_sort);

    // optimize coati FST
    VectorFstStdArc evo_fst, tmp;
    tmp = VectorFstStdArc(evo_comp);
    evo_fst = optimize(tmp);

    VectorFstStdArc evo_rmep;
    evo_rmep = fst::RmEpsilonFst<StdArc>(evo_fst);  // epsilon removal

    return evo_rmep;
}

/// @private
// GCOVR_EXCL_START
TEST_CASE("fst_alignment") {
    coati::alignment_t aln;
    aln.data.path = "test-fst.fasta";

    std::ofstream out;
    out.open("test-fst.fasta");
    REQUIRE(out);
    out << ">1\nCTCTGGATAGTG\n>2\nCTATAGTG\n";
    out.close();

    SUBCASE("coati model, output fasta") {
        aln.model = "coati";
        aln.weight_file = "score.log";
        aln.output = "test-fst-alignment.fasta";

        REQUIRE(fst_alignment(aln));
        std::ifstream infile(aln.output);  // input file
        std::string s1, s2;

        infile >> s1 >> s2;
        CHECK_EQ(s1, ">1");
        CHECK_EQ(s2, "CTCTGGATAGTG");

        infile >> s1 >> s2;
        CHECK_EQ(s1, ">2");
        CHECK_EQ(s2, "CT----ATAGTG");

        std::ifstream inweight(aln.weight_file);
        std::string s;
        inweight >> s;
        REQUIRE(std::filesystem::remove(aln.weight_file));
        REQUIRE(std::filesystem::remove(aln.output));
        CHECK_EQ(s, "test-fst.fasta,coati,9.31397");
    }

    SUBCASE("coati model, output phylip") {
        aln.model = "coati";
        aln.output = "test-fst-phylip.phy";

        REQUIRE(fst_alignment(aln));

        std::ifstream infile(aln.output);
        std::string s1, s2;

        infile >> s1 >> s2;
        CHECK_EQ(s1, "2");
        CHECK_EQ(s2, "12");

        infile >> s1 >> s2;
        CHECK_EQ(s1, "1");
        CHECK_EQ(s2, "CTCTGGATAGTG");

        infile >> s1 >> s2;
        CHECK_EQ(s1, "2");
        CHECK_EQ(s2, "CT----ATAGTG");

        REQUIRE(std::filesystem::remove(aln.output));
    }

    SUBCASE("dna model") {
        aln.model = "dna";
        aln.weight_file = "score.log";
        aln.output = "test-fst-alignment.fasta";

        REQUIRE(fst_alignment(aln));

        std::ifstream infile(aln.output);
        std::string s1, s2;

        infile >> s1 >> s2;
        CHECK_EQ(s1, ">1");
        CHECK_EQ(s2, "CTCTGGATAGTG");

        infile >> s1 >> s2;
        CHECK_EQ(s1, ">2");
        CHECK_EQ(s2, "CT----ATAGTG");
        REQUIRE(std::filesystem::remove(aln.output));

        std::ifstream inweight(aln.weight_file);
        std::string s;
        inweight >> s;
        REQUIRE(std::filesystem::remove(aln.weight_file));
        CHECK_EQ(s, "test-fst.fasta,dna,9.31994");
    }

    SUBCASE("ecm model") {
        aln.model = "ecm";
        aln.weight_file = "score.log";
        aln.output = "test-fst-alignment-ecm.fasta";

        REQUIRE(fst_alignment(aln));
        std::ifstream infile(aln.output);
        std::string s1, s2;

        infile >> s1 >> s2;
        CHECK_EQ(s1, ">1");
        CHECK_EQ(s2, "CTCTGGATAGTG");

        infile >> s1 >> s2;
        CHECK_EQ(s1, ">2");
        CHECK_EQ(s2, "CT----ATAGTG");

        REQUIRE(std::filesystem::remove(aln.output));

        std::ifstream inweight(aln.weight_file);
        std::string s;
        inweight >> s;
        REQUIRE(std::filesystem::remove(aln.weight_file));
        CHECK_EQ(s, "test-fst.fasta,ecm,9.31388");
    }

    SUBCASE("Unknown model") {
        aln.model = "unknown";
        aln.weight_file = "score.log";
        aln.output = "test-fst-alignment.fasta";

        CHECK_THROWS_AS(coati::utils::set_subst(aln), std::invalid_argument);
    }
    SUBCASE("Three sequences - fail") {
        std::ofstream out;
        out.open("test-fst.fasta");
        REQUIRE(out);
        out << ">1\nCTCTGGATAGTG\n>2\nCTATAGTG\n>3\nCTATAGTGTG\n";
        out.close();

        CHECK_THROWS_AS(coati::fst_alignment(aln), std::invalid_argument);
    }
    SUBCASE("Scoring - fail") {
        aln.score = true;
        CHECK_THROWS_AS(coati::fst_alignment(aln), std::invalid_argument);
    }
    SUBCASE("Length of ref not multiple of 3 - fail") {
        std::ofstream out;
        out.open("test-fst.fasta");
        REQUIRE(out);
        out << ">1\nCTCTGGATAGT\n>2\nCTATAGTG\n";
        out.close();

        CHECK_THROWS_AS(coati::fst_alignment(aln), std::invalid_argument);
    }

    REQUIRE(std::filesystem::remove("test-fst.fasta"));
}
// GCOVR_EXCL_STOP

}  // namespace coati
