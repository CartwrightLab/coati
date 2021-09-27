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

#include <coati/align.hpp>
#include <filesystem>

/* Alignment using dynamic programming implementation of marginal COATi model */
bool mcoati(coati::utils::args_t& args, coati::utils::alignment_t& aln) {
    coati::Matrixf P(64, 64), p_marg;
    std::ofstream out_w;

    // score alignment
    if(args.score) {
        std::cout << alignment_score(args, aln.subst_matrix) << std::endl;
        return true;
    }

    size_t len_a = args.fasta.seqs[0].length();
    if((len_a % 3 != 0) && (len_a % args.gap.len != 0)) {
        throw std::invalid_argument(
            "Length of reference sequence must be multiple of 3.");
    }
    if(args.fasta.seqs[1].length() % args.gap.len != 0) {
        throw std::invalid_argument(
            "Length of descendant sequence must be multiple of " +
            std::to_string(args.gap.len) + ".");
    }

    auto anc = args.fasta.seqs[0];
    auto des = args.fasta.seqs[1];
    coati::align_pair_work_t work;
    sequence_pair_t seq_pair = coati::utils::marginal_seq_encoding(anc, des);
    coati::align_pair(work, seq_pair[0], seq_pair[1], aln.subst_matrix, args);
    coati::traceback(work, anc, des, aln, args.gap.len);

    if(!args.weight_file.empty()) {  // save weight and filename
        out_w.open(args.weight_file, std::ios::app | std::ios::out);
        out_w << args.fasta.path << "," << args.model << "," << aln.weight
              << std::endl;
        out_w.close();
    }

    // write alignment
    if(aln.fasta.path.extension() == ".fasta" ||
       aln.fasta.path.extension() == ".fa") {
        return coati::write_fasta(aln.fasta);
    }
    return coati::write_phylip(aln.fasta);
}

TEST_CASE("mcoati") {
    SUBCASE("Alignment - output fasta") {
        coati::utils::args_t args("", {"1", "2"}, {"CTCTGGATAGTG", "CTATAGTG"},
                                  "m-coati", "score.log",
                                  "test-mcoati-fasta.fasta");
        coati::fasta_t result(args.output);
        coati::utils::alignment_t aln;
        aln.fasta.path = args.output;
        aln.fasta.names = args.fasta.names;
        coati::utils::set_subst(args, aln);

        if(std::filesystem::exists(args.output)) {
            std::filesystem::remove(args.output);
        }
        if(std::filesystem::exists(args.weight_file)) {
            std::filesystem::remove(args.weight_file);
        }

        REQUIRE(mcoati(args, aln));
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
        CHECK(s.substr(s.length() - 7) == "1.51294");
    }

    SUBCASE("Alignment - output phylip") {
        coati::utils::args_t args("", {"1", "2"}, {"GCGACTGTT", "GCGATTGCTGTT"},
                                  "m-coati", "", "test-mcoati-phylip.phy");
        coati::fasta_t result(args.output);
        coati::utils::alignment_t aln;
        aln.fasta.path = args.output;
        aln.fasta.names = args.fasta.names;
        coati::utils::set_subst(args, aln);

        if(std::filesystem::exists(args.output)) {
            std::filesystem::remove(args.output);
        }

        REQUIRE(mcoati(args, aln));

        std::ifstream infile(args.output);
        std::string s1, s2;

        infile >> s1 >> s2;
        CHECK(s1.compare("2") == 0);
        CHECK(s2.compare("12") == 0);

        infile >> s1 >> s2;
        CHECK(s1.compare("1") == 0);
        CHECK(s2.compare("GCGA---CTGTT") == 0);

        infile >> s1 >> s2;
        CHECK(s1.compare("2") == 0);
        CHECK(s2.compare("GCGATTGCTGTT") == 0);

        CHECK(std::filesystem::remove(args.output));
    }

    SUBCASE("Alignment 2 dels - output phylip") {
        coati::utils::args_t args("", {"1", "2"}, {"ACGTTAAGGGGT", "ACGAAT"},
                                  "m-coati", "", "test-mcoati-phylip2.phy");
        coati::fasta_t result(args.output);
        coati::utils::alignment_t aln;
        aln.fasta.path = args.output;
        aln.fasta.names = args.fasta.names;
        coati::utils::set_subst(args, aln);

        if(std::filesystem::exists(args.output)) {
            std::filesystem::remove(args.output);
        }

        REQUIRE(mcoati(args, aln));

        std::ifstream infile(args.output);
        std::string s1, s2;

        infile >> s1 >> s2;
        CHECK(s1.compare("2") == 0);
        CHECK(s2.compare("12") == 0);

        infile >> s1 >> s2;
        CHECK(s1.compare("1") == 0);
        CHECK(s2.compare("ACGTTAAGGGGT") == 0);

        infile >> s1 >> s2;
        CHECK(s1.compare("2") == 0);
        CHECK(s2.compare("ACG--AA----T") == 0);

        CHECK(std::filesystem::remove(args.output));
    }

    SUBCASE("Alignment with gap length multiple of 3") {
        coati::utils::args_t args(
            "", {"1", "2"}, {"ACGTTAAGGGGT", "ACGAAT"}, "m-coati", "",
            "test-mcoati-no-frameshifts.fasta", false, "", "", "", 3);
        coati::fasta_t result(args.output);
        coati::utils::alignment_t aln;
        aln.fasta.path = args.output;
        aln.fasta.names = args.fasta.names;
        coati::utils::set_subst(args, aln);

        if(std::filesystem::exists(args.output)) {
            std::filesystem::remove(args.output);
        }

        REQUIRE(mcoati(args, aln));
        result = coati::read_fasta(args.output.string());

        CHECK(std::filesystem::remove(args.output));

        CHECK(result.names[0] == "1");
        CHECK(result.names[1] == "2");

        CHECK(result.seqs[0] == "ACG---TTAAGGGGT");
        CHECK(result.seqs[1] == "ACGAAT---------");
    }

    SUBCASE("Alignment with gap length multiple of 3 - fail") {
        coati::utils::args_t args(
            "", {"1", "2"}, {"GCGATTGCTGT", "GCGACTGTT"}, "m-coati", "",
            "test-mcoati-no-frameshifts.fasta", false, "", "", "", 3);
        coati::utils::alignment_t aln;
        aln.fasta.path = args.output;
        aln.fasta.names = args.fasta.names;
        coati::utils::set_subst(args, aln);
        REQUIRE_THROWS_AS(mcoati(args, aln), std::invalid_argument);
    }

    SUBCASE("Score alignment") {
        coati::utils::args_t args("", {"1", "2"},
                                  {"CTCTGGATAGTG", "CT----ATAGTG"}, "m-coati",
                                  "", "test-mcoati-score.fasta", true);
        coati::utils::alignment_t aln;
        aln.fasta.path = args.output;
        aln.fasta.names = args.fasta.names;
        coati::utils::set_subst(args, aln);
        REQUIRE(mcoati(args, aln));
    }

    SUBCASE("Score alignment - fail") {
        coati::utils::args_t args("", {"1", "2"}, {"CTCTGGATAGTG", "CTATAGTG"},
                                  "m-coati", "", "test-mcoati-score.fasta",
                                  true);
        coati::fasta_t result(args.output);
        coati::utils::alignment_t aln;
        aln.fasta.path = args.output;
        aln.fasta.names = args.fasta.names;
        coati::utils::set_subst(args, aln);

        REQUIRE_THROWS_AS(mcoati(args, aln), std::invalid_argument);
    }
}

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
    if(aln.fasta.path.extension() == ".fasta") {
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

/* Initial msa by collapsing indels after pairwise aln with reference */
bool ref_indel_alignment(coati::utils::args_t& args) {
    coati::Matrixf P(64, 64), p_marg;
    tree_t tree;
    std::string newick;
    coati::utils::alignment_t aln, aln_tmp;

    aln.fasta = coati::fasta_t(args.output);

    // read newick tree file
    if(!read_newick(args.tree, newick)) {
        throw std::invalid_argument("Reading newick tree failed.");
    }

    // parse tree into tree_t (vector<node_t>) variable
    if(parse_newick(newick, tree) != 0) {
        throw std::invalid_argument("Parsing newick tree failed.");
    }

    // reroot tree
    if(!reroot(tree, args.ref)) {
        throw std::invalid_argument("Re-rooting tree failed.");
    }

    // find position of ref in tree
    std::size_t ref_pos = 0;
    if(!find_node(tree, args.ref, ref_pos)) {
        throw std::invalid_argument("Reference node not found in tree.");
    }

    // find sequence of ref in args
    std::vector<std::string> pair_seqs;
    std::string ref_seq;
    if(!find_seq(args.ref, args.fasta, ref_seq)) {
        throw std::invalid_argument("reference sequence " + args.ref +
                                    " not found in fasta file.");
    }

    pair_seqs.push_back(ref_seq);
    pair_seqs.push_back(ref_seq);

    // vector to store insertion_data_t for each node in tree
    std::vector<insertion_data_t> nodes_ins(tree.size());

    // add insertion_data for REF
    nodes_ins[ref_pos] = insertion_data_t(
        ref_seq, args.ref,
        SparseVectorInt(static_cast<Eigen::Index>(2 * ref_seq.length())));

    // pairwise alignment for each leaf
    std::string node_seq;
    for(std::size_t node = 0; node < tree.size(); node++) {
        if(tree[node].is_leaf && (tree[node].label != args.ref)) {
            float branch = distance_ref(tree, ref_pos, node);
            if(!find_seq(tree[node].label, args.fasta, node_seq)) {
                throw std::invalid_argument("sequence " + tree[node].label +
                                            " not found in fasta file.");
            }

            pair_seqs[1] = node_seq;

            args.br_len = branch;
            coati::utils::set_subst(args, aln_tmp);

            aln_tmp.fasta.seqs.clear();
            auto anc = pair_seqs[0];
            auto des = pair_seqs[1];
            coati::align_pair_work_t work;
            sequence_pair_t seq_pair =
                coati::utils::marginal_seq_encoding(anc, des);
            coati::align_pair(work, seq_pair[0], seq_pair[1],
                              aln_tmp.subst_matrix, args);
            coati::traceback(work, anc, des, aln_tmp, args.gap.len);

            SparseVectorInt ins_vector(
                static_cast<Eigen::Index>(aln_tmp.fasta.seqs[1].length()));
            insertion_flags(aln_tmp.fasta.seqs[0], aln_tmp.fasta.seqs[1],
                            ins_vector);

            nodes_ins[node] = insertion_data_t(aln_tmp.fasta.seqs[1],
                                               tree[node].label, ins_vector);
        }
    }

    // get position of inodes in tree and set leafs as visited (true)
    std::vector<std::size_t> inode_indexes;
    std::vector<int> visited(tree.size(), false);  // list of visited nodes

    for(std::size_t node = 0; node < tree.size(); node++) {
        if(!tree[node].is_leaf) {
            inode_indexes.push_back(node);  // add inode position to vector
        } else {
            visited[node] = true;  // set leafs to visited
        }
    }

    // fill list of children
    for(std::size_t i = 0; i < tree.size(); i++) {
        if(tree[i].parent != i) tree[tree[i].parent].children.push_back(i);
    }

    // while not all nodes have been visited (any value in visitied is
    // false)
    while(any_of(visited.begin(), visited.end(), [](bool b) { return !b; })) {
        for(auto inode_pos : inode_indexes) {  // for all inodes
            bool children_visited = true;
            for(auto child : tree[inode_pos].children) {
                if(!visited[child]) {
                    children_visited = false;
                    continue;
                }
            }

            if(!children_visited) {
                continue;  // if all childen of inode have been visited
            }

            visited[inode_pos] = true;

            // if inode only has a child pass information up
            if(tree[inode_pos].children.size() == 1) {
                nodes_ins[inode_pos] = nodes_ins[tree[inode_pos].children[0]];
                continue;
            }

            // create vector of insertion_data_t with children
            std::vector<insertion_data_t> tmp_ins_data(
                tree[inode_pos].children.size());
            for(std::size_t i = 0; i < tree[inode_pos].children.size(); i++) {
                tmp_ins_data[i] = nodes_ins[tree[inode_pos].children[i]];
            }

            // run merge_indels(children_ins_data, nodes_ins[inode_pos]);
            nodes_ins[inode_pos] = insertion_data_t();
            merge_indels(tmp_ins_data, nodes_ins[inode_pos]);
        }
    }

    // transfer result data nodes_ins[ROOT] --> aln && order sequences
    auto root = tree[ref_pos].parent;
    for(const auto& name : args.fasta.names) {
        auto it = find(nodes_ins[root].names.begin(),
                       nodes_ins[root].names.end(), name);
        auto index = distance(nodes_ins[root].names.begin(), it);
        aln.fasta.names.push_back(nodes_ins[root].names[index]);
        aln.fasta.seqs.push_back(nodes_ins[root].sequences[index]);
    }

    // write alignment
    if(std::filesystem::path(aln.fasta.path).extension() == ".fasta") {
        return write_fasta(aln.fasta);
    }
    return coati::write_phylip(aln.fasta);
}

TEST_CASE("ref_indel_alignment") {
    std::ofstream outfile;
    outfile.open("tree-msa.newick");
    REQUIRE(outfile);
    outfile << "((((A:0.1,B:0.1):0.1,C:0.1):0.1,D:0.1):0.1,E:0.1);";
    outfile.close();

    SUBCASE("m-coati model") {
        coati::utils::args_t args(
            "", {"A", "B", "C", "D", "E"},
            {"TCATCG", "TCAGTCG", "TATCG", "TCACTCG", "TCATC"}, "m-coati", "",
            "test-mcoati-msa.fasta", false, "tree-msa.newick", "A");
        coati::fasta_t result(args.output);

        REQUIRE(ref_indel_alignment(args));
        result = coati::read_fasta(args.output.string());

        CHECK(std::filesystem::remove(args.output));

        CHECK(result.names[0] == "A");
        CHECK(result.names[1] == "B");
        CHECK(result.names[2] == "C");
        CHECK(result.names[3] == "D");
        CHECK(result.names[4] == "E");

        CHECK(result.seqs[0] == "TCA--TCG");
        CHECK(result.seqs[1] == "TCA-GTCG");
        CHECK(result.seqs[2] == "T-A--TCG");
        CHECK(result.seqs[3] == "TCAC-TCG");
        CHECK(result.seqs[4] == "TCA--TC-");
    }

    SUBCASE("m-ecm model") {
        coati::utils::args_t args(
            "", {"A", "B", "C", "D", "E"},
            {"TCATCG", "TCAGTCG", "TATCG", "TCACTCG", "TCATC"}, "m-ecm", "",
            "test-mecm-msa.fasta", false, "tree-msa.newick", "A");
        coati::fasta_t result(args.output);

        REQUIRE(ref_indel_alignment(args));
        result = coati::read_fasta(args.output.string());

        CHECK(std::filesystem::remove(args.output));

        CHECK(result.names[0] == "A");
        CHECK(result.names[1] == "B");
        CHECK(result.names[2] == "C");
        CHECK(result.names[3] == "D");
        CHECK(result.names[4] == "E");

        CHECK(result.seqs[0] == "TCA--TCG");
        CHECK(result.seqs[1] == "TCA-GTCG");
        CHECK(result.seqs[2] == "T-A--TCG");
        CHECK(result.seqs[3] == "TCAC-TCG");
        CHECK(result.seqs[4] == "TCA--TC-");
    }
}

float alignment_score(coati::utils::args_t& args, coati::Matrixf& p_marg) {
    std::vector<std::string> aln = args.fasta.seqs;
    if(aln[0].length() != aln[1].length()) {
        throw std::invalid_argument(
            "For alignment scoring both sequences must have equal length.");
    }
    std::string anc{aln[0]};
    boost::erase_all(anc, "-");
    sequence_pair_t seq_pair = coati::utils::marginal_seq_encoding(anc, aln[1]);

    // calculate log(1-g) log(1-e) log(g) log(e) log(pi)
    float_t no_gap = std::log1pf(-args.gap.open);
    float_t gap_stop = std::log1pf(-args.gap.extend);
    float_t gap_open = ::logf(args.gap.open);
    float_t gap_extend = ::logf(args.gap.extend);
    std::vector<coati::float_t> pi{args.pi};
    for(size_t i = 0; i < 4; i++) {
        pi[i] = ::logf(pi[i]);
    }

    float weight{0.f};
    int state{0}, ngap{0};
    for(size_t i = 0; i < aln[0].length(); i++) {
        switch(state) {
        case 0:
            if(aln[0][i] == '-') {
                // insertion;
                weight += gap_open;
                state = 2;
                ngap++;
            } else if(aln[1][i] == '-') {
                // deletion;
                weight += no_gap + gap_open;
                state = 1;
            } else {
                // match/mismatch;
                weight +=
                    2 * no_gap + p_marg(seq_pair[0][i - ngap], seq_pair[1][i]);
            }
            break;
        case 1:
            if(aln[0][i] == '-') {
                throw std::runtime_error(
                    "Insertion after deletion is not modeled.");
            } else if(aln[1][i] == '-') {
                // deletion_ext
                weight += gap_extend;
            } else {
                // match/mismatch
                weight +=
                    gap_stop + p_marg(seq_pair[0][i - ngap], seq_pair[1][i]);
                state = 0;
            }
            break;
        case 2:
            if(aln[0][i] == '-') {
                // insertion_ext
                weight += gap_extend;
                ngap++;
            } else if(aln[1][i] == '-') {
                // deletion
                weight += gap_stop + gap_open;
                state = 1;
            } else {
                // match/mismatch
                weight += gap_stop + no_gap +
                          p_marg(seq_pair[0][i - ngap], seq_pair[1][i]);
                state = 0;
            }
            break;
        }
    }
    // terminal state weight
    if(state == 0) {
        weight += no_gap;
    } else if(state == 2) {
        weight += gap_stop;
    }
    return weight;
}

TEST_CASE("alignment_score") {
    coati::utils::args_t args;
    coati::Matrixf P(mg94_p(0.0133, 0.2, {0.308, 0.185, 0.199, 0.308}));
    coati::Matrixf p_marg = marginal_p(P, args.pi);

    args.fasta.seqs = {"CTCTGGATAGTG", "CT----ATAGTG"};
    REQUIRE(alignment_score(args, p_marg) == doctest::Approx(1.51294f));

    args.fasta.seqs = {"CTCT--AT", "CTCTGGAT"};
    REQUIRE(alignment_score(args, p_marg) == doctest::Approx(-0.835939f));

    args.fasta.seqs = {"CTC", "CT"};
    REQUIRE_THROWS_AS(alignment_score(args, p_marg), std::invalid_argument);
}
