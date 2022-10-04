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
 * \brief Pairwise alignment using FST composition.
 *
 * Composition of indel and subsitution FSTs together with input sequences
 *  as FSAs. Alignment is found searching for shortest path in resulting
 * FST.
 *
 * @param[in] aln coati::alignment_t alignment information.
 */
bool fst_alignment(coati::alignment_t& aln) {
    using fst::StdArc;

    // set substitution matrix according to model
    coati::utils::set_subst(aln);

    // get indel FST
    VectorFstStdArc indel_fst = indel(aln.gap.open, aln.gap.extend, aln.pi);

    // sort mutation and indel FSTs
    VectorFstStdArc mutation_sort, indel_sort;
    mutation_sort = fst::ArcSortFst<StdArc, fst::OLabelCompare<StdArc>>(
        aln.subst_fst, fst::OLabelCompare<StdArc>());
    indel_sort = fst::ArcSortFst<StdArc, fst::ILabelCompare<StdArc>>(
        indel_fst, fst::ILabelCompare<StdArc>());

    // compose mutation and indel FSTs
    fst::ComposeFst<StdArc> coati_comp =
        fst::ComposeFst<StdArc>(mutation_sort, indel_sort);

    // optimize coati FST
    VectorFstStdArc coati_fst;
    coati_fst = optimize(VectorFstStdArc(coati_comp));

    VectorFstStdArc coati_rmep;
    coati_rmep = fst::RmEpsilonFst<StdArc>(coati_fst);  // epsilon removal

    // find alignment graph
    // 1. compose in_tape and coati FSTs
    fst::ComposeFst<StdArc> aln_inter =
        fst::ComposeFst<StdArc>(aln.data.fsts[0], coati_rmep);
    // 2. sort intermediate composition
    VectorFstStdArc aln_inter_sort;
    aln_inter_sort = fst::ArcSortFst<StdArc, fst::OLabelCompare<StdArc>>(
        aln_inter, fst::OLabelCompare<StdArc>());
    // 3. compose intermediate and out_tape FSTs
    VectorFstStdArc graph_fst;
    fst::Compose(aln_inter_sort, aln.data.fsts[1], &graph_fst);

    //  case 1: ComposeFst is time: O(v1 v2 d1 (log d2 + m2)), space O(v1 v2)
    //			ShortestPath with ComposeFst (PDT) is time: O((V+E)^4),
    //          space: O((V+E)^3)
    //          Then convert to VectorFst (FST) is time, space: O(e^(O(V+E)))
    //  case 2: ComposeFst is time: O(v1 v2 d1 (log d2 + m2)), space: O(v1 v2)
    //			Convert ComposeFst (PDT) to VectorFst(FST) is time,
    //          space: O(e^(O(V+E)))
    //          Then ShortestPath with VectorFst is O(V log(V) + E)
    //  case 3: Compose is time: O(V1 V2 D1 (log D2 + M2)), space O(V1 V2 D1 M2)
    //			Then ShortestPath with VectorFst is O(V log(V) + E)

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
        out_w << aln.data.path << "," << aln.model << "," << distance[0]
              << std::endl;
        out_w.close();
    }

    // topsort path FST
    fst::TopSort(&aln_path);

    // write alignment
    coati::io::write_output(aln.data, aln_path);
    return true;
}

/// @private
// GCOVR_EXCL_START
TEST_CASE("fst_alignment") {
    std::vector<VectorFstStdArc> fsts;
    VectorFstStdArc fsa0, fsa1;

    CHECK(acceptor("CTCTGGATAGTG", fsa0));
    fsts.push_back(fsa0);
    CHECK(acceptor("CTATAGTG", fsa1));
    fsts.push_back(fsa1);

    SUBCASE("coati model, output fasta") {
        coati::alignment_t aln;
        aln.data = coati::data_t("", {"1", "2"}, {"CTCTGGATAGTG", "CTATAGTG"});
        aln.model = "coati";
        aln.weight_file = "score.log";
        aln.data.out_file = {{"test-fst-alignment.fasta"}, {".fasta"}};
        aln.data.fsts = fsts;

        if(std::filesystem::exists(aln.data.out_file.path)) {
            std::filesystem::remove(aln.data.out_file.path);
        }
        if(std::filesystem::exists(aln.weight_file)) {
            std::filesystem::remove(aln.weight_file);
        }

        REQUIRE(fst_alignment(aln));
        std::ifstream infile(aln.data.out_file.path);  // input file
        std::string s1, s2;

        infile >> s1 >> s2;
        CHECK_EQ(s1, ">1");
        CHECK_EQ(s2, "CTCTGGATAGTG");

        infile >> s1 >> s2;
        CHECK_EQ(s1, ">2");
        CHECK_EQ(s2, "CT----ATAGTG");
        // CHECK(std::filesystem::remove(aln.data.out_file.path));

        std::ifstream inweight(aln.weight_file);
        std::string s;
        inweight >> s;
        CHECK(std::filesystem::remove(aln.weight_file));
        CHECK_EQ(s.substr(s.length() - 7), "9.31397");
    }

    SUBCASE("coati model, output phylip") {
        coati::alignment_t aln;
        aln.data = coati::data_t("", {"1", "2"}, {"CTCTGGATAGTG", "CTATAGTG"});
        aln.model = "coati";
        aln.weight_file = "";
        aln.data.out_file = {{"test-fst-phylip.phy"}, {".phy"}};
        aln.data.fsts = fsts;

        if(std::filesystem::exists(aln.data.out_file.path)) {
            std::filesystem::remove(aln.data.out_file.path);
        }
        if(std::filesystem::exists(aln.weight_file)) {
            std::filesystem::remove(aln.weight_file);
        }

        REQUIRE(fst_alignment(aln));

        std::ifstream infile(aln.data.out_file.path);
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

        CHECK(std::filesystem::remove(aln.data.out_file.path));
    }

    SUBCASE("dna model") {
        coati::alignment_t aln;
        aln.data = coati::data_t("", {"1", "2"}, {"CTCTGGATAGTG", "CTATAGTG"});
        aln.model = "dna";
        aln.weight_file = "score.log";
        aln.data.out_file = {{"test-fst-alignment.fasta"}, {".fasta"}};
        aln.data.fsts = fsts;

        if(std::filesystem::exists(aln.data.out_file.path)) {
            std::filesystem::remove(aln.data.out_file.path);
        }
        if(std::filesystem::exists(aln.weight_file)) {
            std::filesystem::remove(aln.weight_file);
        }

        REQUIRE(fst_alignment(aln));

        std::ifstream infile(aln.data.out_file.path);
        std::string s1, s2;

        infile >> s1 >> s2;
        CHECK_EQ(s1, ">1");
        CHECK_EQ(s2, "CTCTGGATAGTG");

        infile >> s1 >> s2;
        CHECK_EQ(s1, ">2");
        CHECK_EQ(s2, "CT----ATAGTG");
        CHECK(std::filesystem::remove(aln.data.out_file.path));

        std::ifstream inweight(aln.weight_file);
        std::string s;
        inweight >> s;
        CHECK(std::filesystem::remove(aln.weight_file));
        CHECK_EQ(s.substr(s.length() - 7), "9.31994");
    }

    SUBCASE("ecm model") {
        coati::alignment_t aln;
        aln.data = coati::data_t("", {"1", "2"}, {"CTCTGGATAGTG", "CTATAGTG"});
        aln.model = "ecm";
        aln.weight_file = "score.log";
        aln.data.out_file = {{"test-fst-alignment-ecm.fasta"}, {".fasta"}};
        aln.data.fsts = fsts;

        if(std::filesystem::exists(aln.data.out_file.path)) {
            std::filesystem::remove(aln.data.out_file.path);
        }
        if(std::filesystem::exists(aln.weight_file)) {
            std::filesystem::remove(aln.weight_file);
        }

        REQUIRE(fst_alignment(aln));
        std::ifstream infile(aln.data.out_file.path);
        std::string s1, s2;

        infile >> s1 >> s2;
        CHECK_EQ(s1, ">1");
        CHECK_EQ(s2, "CTCTGGATAGTG");

        infile >> s1 >> s2;
        CHECK_EQ(s1, ">2");
        CHECK_EQ(s2, "CT----ATAGTG");

        CHECK(std::filesystem::remove(aln.data.out_file.path));

        std::ifstream inweight(aln.weight_file);
        std::string s;
        inweight >> s;
        CHECK(std::filesystem::remove(aln.weight_file));
        CHECK_EQ(s.substr(s.length() - 7), "9.31388");
    }

    SUBCASE("Unknown model") {
        coati::alignment_t aln;
        aln.data = coati::data_t("", {"1", "2"}, {"CTCTGGATAGTG", "CTATAGTG"});
        aln.model = "unknown";
        aln.weight_file = "score.log";
        aln.data.out_file = {{"test-fst-alignment.fasta"}, {".fasta"}};

        REQUIRE_THROWS_AS(coati::utils::set_subst(aln), std::invalid_argument);
    }
}
// GCOVR_EXCL_STOP

}  // namespace coati
