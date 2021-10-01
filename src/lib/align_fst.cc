/*
# Copyright (c) 2020-2021 Juan J. Garcia Mesa <juanjosegarciamesa@gmail.com>
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

/* Alignment using FST library*/
bool fst_alignment(coati::utils::args_t& args, coati::utils::alignment_t& aln) {
    using fst::StdArc;

    // get indel FST
    VectorFstStdArc indel_fst = indel(args.gap.open, args.gap.extend, args.pi);

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
        fst::ComposeFst<StdArc>(aln.seqs[0], coati_rmep);
    // 2. sort intermediate composition
    VectorFstStdArc aln_inter_sort;
    aln_inter_sort = fst::ArcSortFst<StdArc, fst::OLabelCompare<StdArc>>(
        aln_inter, fst::OLabelCompare<StdArc>());
    // 3. compose intermediate and out_tape FSTs
    VectorFstStdArc graph_fst;
    fst::Compose(aln_inter_sort, aln.seqs[1], &graph_fst);

    //  case 1: ComposeFst is time: O(v1 v2 d1 (log d2 + m2)), space O(v1 v2)
    //			ShortestPath with ComposeFst (PDT) is time: O((V+E)^3),
    //          space: O((V+E)^3)
    //          Then convert to VectorFst (FST) is time, space: O(e^(O(V+E)))
    //  case 2: ComposeFst is time: O(v1 v2 d1 (log d2 + m2)), space: O(v1 v2)
    //			Convert ComposeFst (PDT) to VectorFst(FST) is time,
    //          space: O(e^(O(V+E)))
    // 			Then ShortestPath with VectorFst is O(V log(V) + E)
    //  case 3: Compose is time: O(V1 V2 D1 (log D2 + M2)), space O(V1 V2 D1 M2)
    //			Then ShortestPath with VectorFst is O(V log(V) + E)

    // find shortest path through graph
    VectorFstStdArc aln_path;
    fst::ShortestPath(graph_fst, &aln_path);

    // shortestdistance = weight of shortestpath
    if(!args.weight_file.empty()) {
        std::vector<StdArc::Weight> distance;
        std::ofstream out_w;

        fst::ShortestDistance(aln_path, &distance);
        // append weight and fasta file name info in file
        out_w.open(args.weight_file, std::ios::app | std::ios::out);
        out_w << args.fasta.path << "," << args.model << "," << distance[0]
              << std::endl;
        out_w.close();
    }

    // topsort path FST
    fst::TopSort(&aln_path);

    // write alignment
    if((aln.fasta.path.extension() == ".fasta") ||
       (aln.fasta.path.extension() == ".fa")) {
        return coati::write_fasta(aln_path, aln.fasta);
    }
    return coati::write_phylip(aln_path, aln.fasta);
}

TEST_CASE("fst_alignment") {
    std::vector<VectorFstStdArc> fsts;
    VectorFstStdArc fsa0, fsa1;

    CHECK(acceptor("CTCTGGATAGTG", fsa0));
    fsts.push_back(fsa0);
    CHECK(acceptor("CTATAGTG", fsa1));
    fsts.push_back(fsa1);

    SUBCASE("coati model, output fasta") {
        coati::utils::args_t args("", {"1", "2"}, {"CTCTGGATAGTG", "CTATAGTG"},
                                  "coati", "score.log",
                                  "test-fst-alignment.fasta");
        coati::fasta_t result(args.output);
        coati::utils::alignment_t aln;
        aln.fasta.path = args.output;
        aln.fasta.names = args.fasta.names;
        coati::utils::set_subst(args, aln);
        aln.seqs = fsts;

        if(std::filesystem::exists(args.output)) {
            std::filesystem::remove(args.output);
        }
        if(std::filesystem::exists(args.weight_file)) {
            std::filesystem::remove(args.weight_file);
        }

        REQUIRE(fst_alignment(args, aln));
        result = coati::read_fasta(args.output.string());

        CHECK(std::filesystem::remove(args.output));

        CHECK(result.names[0] == "1");
        CHECK(result.names[1] == "2");

        CHECK(result.seqs[0] == "CTCTGGATAGTG");
        CHECK(result.seqs[1] == "CT----ATAGTG");

        std::ifstream infile(args.weight_file);
        std::string s;
        infile >> s;
        CHECK(std::filesystem::remove(args.weight_file));
        CHECK(s.substr(s.length() - 7) == "9.31397");
    }

    SUBCASE("coati model, output phylip") {
        coati::utils::args_t args("", {"1", "2"}, {"CTCTGGATAGTG", "CTATAGTG"},
                                  "coati", "", "test-fst-phylip.phy");
        coati::fasta_t result(args.output);
        coati::utils::alignment_t aln;
        aln.fasta.path = args.output;
        aln.fasta.names = args.fasta.names;
        coati::utils::set_subst(args, aln);
        aln.seqs = fsts;

        if(std::filesystem::exists(args.output)) {
            std::filesystem::remove(args.output);
        }
        if(std::filesystem::exists(args.weight_file)) {
            std::filesystem::remove(args.weight_file);
        }

        REQUIRE(fst_alignment(args, aln));

        std::ifstream infile(args.output);
        std::string s1, s2;

        infile >> s1 >> s2;
        CHECK(s1.compare("2") == 0);
        CHECK(s2.compare("12") == 0);

        infile >> s1 >> s2;
        CHECK(s1.compare("1") == 0);
        CHECK(s2.compare("CTCTGGATAGTG") == 0);

        infile >> s1 >> s2;
        CHECK(s1.compare("2") == 0);
        CHECK(s2.compare("CT----ATAGTG") == 0);

        CHECK(std::filesystem::remove(args.output));
    }

    SUBCASE("dna model") {
        coati::utils::args_t args("", {"1", "2"}, {"CTCTGGATAGTG", "CTATAGTG"},
                                  "dna", "score.log",
                                  "test-fst-alignment.fasta");
        coati::fasta_t result(args.output);
        coati::utils::alignment_t aln;
        aln.fasta.path = args.output;
        aln.fasta.names = args.fasta.names;
        coati::utils::set_subst(args, aln);
        aln.seqs = fsts;

        if(std::filesystem::exists(args.output)) {
            std::filesystem::remove(args.output);
        }
        if(std::filesystem::exists(args.weight_file)) {
            std::filesystem::remove(args.weight_file);
        }

        REQUIRE(fst_alignment(args, aln));
        result = coati::read_fasta(args.output.string());

        CHECK(std::filesystem::remove(args.output));

        CHECK(result.names[0] == "1");
        CHECK(result.names[1] == "2");

        CHECK(result.seqs[0] == "CTCTGGATAGTG");
        CHECK(result.seqs[1] == "CT----ATAGTG");

        std::ifstream infile(args.weight_file);
        std::string s;
        infile >> s;
        CHECK(std::filesystem::remove(args.weight_file));
        std::cout << s.substr(s.length() - 7) << std::endl;
        CHECK(s.substr(s.length() - 7) == "9.31994");
    }

    SUBCASE("ecm model") {
        coati::utils::args_t args("", {"1", "2"}, {"CTCTGGATAGTG", "CTATAGTG"},
                                  "ecm", "score.log",
                                  "test-fst-alignment.fasta");
        coati::fasta_t result(args.output);
        coati::utils::alignment_t aln;
        aln.fasta.path = args.output;
        aln.fasta.names = args.fasta.names;
        coati::utils::set_subst(args, aln);
        aln.seqs = fsts;

        if(std::filesystem::exists(args.output)) {
            std::filesystem::remove(args.output);
        }
        if(std::filesystem::exists(args.weight_file)) {
            std::filesystem::remove(args.weight_file);
        }

        REQUIRE(fst_alignment(args, aln));
        result = coati::read_fasta(args.output.string());

        CHECK(std::filesystem::remove(args.output));

        CHECK(result.names[0] == "1");
        CHECK(result.names[1] == "2");

        CHECK(result.seqs[0] == "CTCTGGATAGTG");
        CHECK(result.seqs[1] == "CT----ATAGTG");

        std::ifstream infile(args.weight_file);
        std::string s;
        infile >> s;
        CHECK(std::filesystem::remove(args.weight_file));
        CHECK(s.substr(s.length() - 7) == "9.31388");
    }

    SUBCASE("Unknown model") {
        coati::utils::args_t args("", {"1", "2"}, {"CTCTGGATAGTG", "CTATAGTG"},
                                  "unknown", "", "test-fst-alignment.fasta");
        coati::utils::alignment_t aln;
        aln.fasta.path = args.output;
        aln.fasta.names = args.fasta.names;

        REQUIRE_THROWS_AS(coati::utils::set_subst(args, aln),
                          std::invalid_argument);
    }
}

}  // namespace coati