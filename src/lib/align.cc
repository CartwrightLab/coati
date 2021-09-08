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
int mcoati(input_t& in_data) {
    Matrix P(64, 64);
    std::vector<VectorFstStdArc> fsts;
    std::ofstream out_w;
    alignment_t aln;
    aln.f.seq_names = in_data.fasta_file.seq_names;
    aln.f.path = in_data.out_file;

    if(!in_data.rate.empty()) {
        in_data.mut_model = "user_marg_model";
        P = parse_matrix_csv(in_data.rate);
    } else if(in_data.mut_model.compare("m-coati") == 0) {
        P = mg94_p(in_data.br_len, in_data.omega);
    } else {  // m-ecm
        P = ecm_p(in_data.br_len, in_data.omega);
    }

    if(in_data.score) {
        std::cout << alignment_score(in_data.fasta_file.seq_data, P,
                                     in_data.gapo, in_data.gape)
                  << std::endl;
        return EXIT_SUCCESS;
    }

    if(in_data.frameshifts) {
        if(mg94_marginal(in_data.fasta_file.seq_data, aln, P) != 0) {
            return EXIT_FAILURE;
        }
    } else {
        if(mg94_marginal(in_data.fasta_file.seq_data, aln, P, false) != 0) {
            return EXIT_FAILURE;
        }
    }

    if(!in_data.weight_file.empty()) {
        // append weight and fasta file name to file
        out_w.open(in_data.weight_file, std::ios::app | std::ios::out);
        out_w << in_data.fasta_file.path << "," << in_data.mut_model << ","
              << aln.weight << std::endl;
        out_w.close();
    }

    // write alignment
    if(aln.f.path.extension() == ".fasta" || aln.f.path.extension() == ".fa") {
        return write_fasta(aln.f);
    }
    return write_phylip(aln.f);
}

TEST_CASE("mcoati") {
    Matrix P(mg94_p(0.0133, 0.2));

    SUBCASE("Alignment with frameshifts (default) - output fasta") {
        input_t input_data("", {"1", "2"}, {"CTCTGGATAGTG", "CTATAGTG"},
                           "m-coati", "score.log", "test-mcoati-fasta.fasta");
        fasta_t result(input_data.out_file);

        if(std::filesystem::exists(input_data.out_file)) {
            std::filesystem::remove(input_data.out_file);
        }
        if(std::filesystem::exists(input_data.weight_file)) {
            std::filesystem::remove(input_data.weight_file);
        }

        REQUIRE(mcoati(input_data) == 0);
        REQUIRE(read_fasta(result) == 0);

        CHECK(std::filesystem::remove(input_data.out_file));

        CHECK(result.seq_names[0] == "1");
        CHECK(result.seq_names[1] == "2");

        CHECK(result.seq_data[0] == "CTCTGGATAGTG");
        CHECK(result.seq_data[1] == "CT----ATAGTG");

        std::ifstream infile(input_data.weight_file);
        std::string s;
        infile >> s;
        CHECK(std::filesystem::remove(input_data.weight_file));
        CHECK(s.substr(s.length() - 7) == "9.29164");
    }

    SUBCASE("Alignment with frameshifts (default) - output phylip") {
        input_t input_data("", {"1", "2"}, {"GCGACTGTT", "GCGATTGCTGTT"},
                           "m-coati", "score.log", "test-mcoati-phylip.phy");
        fasta_t result(input_data.out_file);

        if(std::filesystem::exists(input_data.out_file)) {
            std::filesystem::remove(input_data.out_file);
        }

        REQUIRE(mcoati(input_data) == 0);

        std::ifstream infile(input_data.out_file);
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

        CHECK(std::filesystem::remove(input_data.out_file));
    }

    SUBCASE("Alignment with frameshifts (default) 2 dels - output phylip") {
        input_t input_data("", {"1", "2"}, {"ACGTTAAGGGGT", "ACGAAT"},
                           "m-coati", "score.log", "test-mcoati-phylip2.phy");
        fasta_t result(input_data.out_file);

        if(std::filesystem::exists(input_data.out_file)) {
            std::filesystem::remove(input_data.out_file);
        }

        REQUIRE(mcoati(input_data) == 0);

        std::ifstream infile(input_data.out_file);
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

        CHECK(std::filesystem::remove(input_data.out_file));
    }

    SUBCASE("Alignment with no frameshifts") {
        input_t input_data(
            "", {"1", "2"}, {"ACGTTAAGGGGT", "ACGAAT"}, "m-coati", "score.log",
            "test-mcoati-no-frameshifts.fasta", "", "", "", false, false);
        fasta_t result(input_data.out_file);

        if(std::filesystem::exists(input_data.out_file)) {
            std::filesystem::remove(input_data.out_file);
        }

        REQUIRE(mcoati(input_data) == 0);
        REQUIRE(read_fasta(result) == 0);

        CHECK(std::filesystem::remove(input_data.out_file));

        CHECK(result.seq_names[0] == "1");
        CHECK(result.seq_names[1] == "2");

        CHECK(result.seq_data[0] == "ACG---TTAAGGGGT");
        CHECK(result.seq_data[1] == "ACGAAT---------");
    }

    SUBCASE("No frameshifts length not multiple of 3 - fail") {
        input_t input_data("", {"1", "2"}, {"GCGATTGCTGT", "GCGACTGTT"},
                           "no_frameshifts", "score.log",
                           "test-mcoati-no-frameshifts.fasta");
        REQUIRE_THROWS_AS(mcoati(input_data), std::invalid_argument);
    }

    SUBCASE("Score alignment") {
        input_t input_data("", {"1", "2"}, {"CTCTGGATAGTG", "CT----ATAGTG"},
                           "m-coati", "score.log", "test-mcoati-score.fasta",
                           "", "", "", true);
        REQUIRE(mcoati(input_data) == 0);
    }

    SUBCASE("Reference length not multiple of 3 - fail") {
        input_t input_data("", {"1", "2"}, {"CTCTGGATAGT", "CT----ATAGTG"},
                           "m-coati", "score.log", "test-mcoati-score.fasta",
                           "", "", "", true);
        REQUIRE_THROWS_AS(mcoati(input_data), std::invalid_argument);
    }

    SUBCASE("Score alignment - fail") {
        input_t input_data("", {"1", "2"}, {"CTCTGGATAGTG", "CTATAGTG"},
                           "m-coati", "score.log", "test-mcoati-score.fasta",
                           "", "", "", true);
        fasta_t result(input_data.out_file);

        REQUIRE_THROWS_AS(mcoati(input_data), std::invalid_argument);
    }
}

/* Alignment using FST library*/
int fst_alignment(input_t& in_data, std::vector<VectorFstStdArc>& fsts) {
    using fst::StdArc;

    VectorFstStdArc mut_fst;

    if(in_data.mut_model.compare("coati") == 0) {
        mut_fst = mg94(in_data.br_len, in_data.omega);
    } else if(in_data.mut_model.compare("dna") == 0) {
        mut_fst = dna(in_data.br_len, in_data.omega);
    } else if(in_data.mut_model.compare("ecm") == 0) {
        mut_fst = ecm(in_data.br_len, in_data.omega);
    } else {
        throw std::invalid_argument("Mutation model unknown. Exiting!");
    }

    // get indel FST
    VectorFstStdArc indel_fst =
        indel(in_data.mut_model, in_data.gapo, in_data.gape);

    // sort mutation and indel FSTs
    VectorFstStdArc mutation_sort, indel_sort;
    mutation_sort = fst::ArcSortFst<StdArc, fst::OLabelCompare<StdArc>>(
        mut_fst, fst::OLabelCompare<StdArc>());
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
        fst::ComposeFst<StdArc>(fsts[0], coati_rmep);
    // 2. sort intermediate composition
    VectorFstStdArc aln_inter_sort;
    aln_inter_sort = fst::ArcSortFst<StdArc, fst::OLabelCompare<StdArc>>(
        aln_inter, fst::OLabelCompare<StdArc>());
    // 3. compose intermediate and out_tape FSTs
    VectorFstStdArc graph_fst;
    fst::Compose(aln_inter_sort, fsts[1], &graph_fst);

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
    if(!in_data.weight_file.empty()) {
        std::vector<StdArc::Weight> distance;
        std::ofstream out_w;

        fst::ShortestDistance(aln_path, &distance);
        // append weight and fasta file name info in file
        out_w.open(in_data.weight_file, std::ios::app | std::ios::out);
        out_w << in_data.fasta_file.path << "," << in_data.mut_model << ","
              << distance[0] << std::endl;
        out_w.close();
    }

    // topsort path FST
    fst::TopSort(&aln_path);

    fasta_t out_fasta(in_data.out_file, in_data.fasta_file.seq_names);

    // write alignment
    if(out_fasta.path.extension() == ".fasta") {
        return write_fasta(aln_path, out_fasta);
    }
    return write_phylip(aln_path, out_fasta);
}

TEST_CASE("fst_alignment") {
    std::vector<VectorFstStdArc> fsts;
    VectorFstStdArc fsa0, fsa1;

    CHECK(acceptor("CTCTGGATAGTG", fsa0));
    fsts.push_back(fsa0);
    CHECK(acceptor("CTATAGTG", fsa1));
    fsts.push_back(fsa1);

    SUBCASE("coati model, output fasta") {
        input_t input_data("", {"1", "2"}, {"CTCTGGATAGTG", "CTATAGTG"},
                           "coati", "score.log", "test-fst-alignment.fasta");
        fasta_t result(input_data.out_file);

        if(std::filesystem::exists(input_data.out_file)) {
            std::filesystem::remove(input_data.out_file);
        }
        if(std::filesystem::exists(input_data.weight_file)) {
            std::filesystem::remove(input_data.weight_file);
        }

        REQUIRE(fst_alignment(input_data, fsts) == 0);
        REQUIRE(read_fasta(result) == 0);

        CHECK(std::filesystem::remove(input_data.out_file));

        CHECK(result.seq_names[0] == "1");
        CHECK(result.seq_names[1] == "2");

        CHECK(result.seq_data[0] == "CTCTGGATAGTG");
        CHECK(result.seq_data[1] == "CT----ATAGTG");

        std::ifstream infile(input_data.weight_file);
        std::string s;
        infile >> s;
        CHECK(std::filesystem::remove(input_data.weight_file));
        CHECK(s.substr(s.length() - 7) == "9.31397");
    }

    SUBCASE("coati model, output phylip") {
        input_t input_data("", {"1", "2"}, {"CTCTGGATAGTG", "CTATAGTG"},
                           "coati", "score.log", "test-fst-phylip.phy");
        fasta_t result(input_data.out_file);

        if(std::filesystem::exists(input_data.out_file)) {
            std::filesystem::remove(input_data.out_file);
        }
        if(std::filesystem::exists(input_data.weight_file)) {
            std::filesystem::remove(input_data.weight_file);
        }

        REQUIRE(fst_alignment(input_data, fsts) == 0);

        std::ifstream infile(input_data.out_file);
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

        CHECK(std::filesystem::remove(input_data.out_file));
        CHECK(std::filesystem::remove(input_data.weight_file));
    }

    SUBCASE("dna model") {
        input_t input_data("", {"1", "2"}, {"CTCTGGATAGTG", "CTATAGTG"}, "dna",
                           "score.log", "test-fst-alignment.fasta");
        fasta_t result(input_data.out_file);

        if(std::filesystem::exists(input_data.out_file)) {
            std::filesystem::remove(input_data.out_file);
        }
        if(std::filesystem::exists(input_data.weight_file)) {
            std::filesystem::remove(input_data.weight_file);
        }

        REQUIRE(fst_alignment(input_data, fsts) == 0);
        REQUIRE(read_fasta(result) == 0);

        CHECK(std::filesystem::remove(input_data.out_file));

        CHECK(result.seq_names[0] == "1");
        CHECK(result.seq_names[1] == "2");

        CHECK(result.seq_data[0] == "CTCTGGATAGTG");
        CHECK(result.seq_data[1] == "CT----ATAGTG");

        std::ifstream infile(input_data.weight_file);
        std::string s;
        infile >> s;
        CHECK(std::filesystem::remove(input_data.weight_file));
        std::cout << s.substr(s.length() - 7) << std::endl;
        CHECK(s.substr(s.length() - 7) == "9.31994");
    }

    SUBCASE("ecm model") {
        input_t input_data("", {"1", "2"}, {"CTCTGGATAGTG", "CTATAGTG"}, "ecm",
                           "score.log", "test-fst-alignment.fasta");
        fasta_t result(input_data.out_file);

        if(std::filesystem::exists(input_data.out_file)) {
            std::filesystem::remove(input_data.out_file);
        }
        if(std::filesystem::exists(input_data.weight_file)) {
            std::filesystem::remove(input_data.weight_file);
        }

        REQUIRE(fst_alignment(input_data, fsts) == 0);
        REQUIRE(read_fasta(result) == 0);

        CHECK(std::filesystem::remove(input_data.out_file));

        CHECK(result.seq_names[0] == "1");
        CHECK(result.seq_names[1] == "2");

        CHECK(result.seq_data[0] == "CTCTGGATAGTG");
        CHECK(result.seq_data[1] == "CT----ATAGTG");

        std::ifstream infile(input_data.weight_file);
        std::string s;
        infile >> s;
        CHECK(std::filesystem::remove(input_data.weight_file));
        CHECK(s.substr(s.length() - 7) == "9.31388");
    }

    SUBCASE("Unknown model") {
        input_t input_data("", {"1", "2"}, {"CTCTGGATAGTG", "CTATAGTG"},
                           "unknown", "score.log", "test-fst-alignment.fasta");

        REQUIRE_THROWS_AS(fst_alignment(input_data, fsts),
                          std::invalid_argument);
    }
}

/* Initial msa by collapsing indels after pairwise aln with reference */
int ref_indel_alignment(input_t& in_data) {
    Matrix P(64, 64);
    tree_t tree;
    std::string newick;
    alignment_t aln, aln_tmp;

    aln.f = fasta_t(in_data.out_file);

    // read newick tree file
    if(!read_newick(in_data.tree, newick)) {
        throw std::invalid_argument("Reading newick tree failed.");
    }

    // parse tree into tree_t (vector<node_t>) variable
    if(parse_newick(newick, tree) != 0) {
        throw std::invalid_argument("Parsing newick tree failed.");
    }

    // reroot tree
    if(!reroot(tree, in_data.ref)) {
        throw std::invalid_argument("Re-rooting tree failed.");
    }

    // find position of ref in tree
    std::size_t ref_pos = 0;
    if(!find_node(tree, in_data.ref, ref_pos)) {
        throw std::invalid_argument("Reference node not found in tree.");
    }

    // find sequence of ref in in_data
    std::vector<std::string> pair_seqs;
    std::string ref_seq;
    if(!find_seq(in_data.ref, in_data.fasta_file, ref_seq)) {
        throw std::invalid_argument("reference sequence " + in_data.ref +
                                    " not found in fasta file.");
    }

    pair_seqs.push_back(ref_seq);
    pair_seqs.push_back(ref_seq);

    // vector to store insertion_data_t for each node in tree
    std::vector<insertion_data_t> nodes_ins(tree.size());

    // add insertion_data for REF
    nodes_ins[ref_pos] = insertion_data_t(
        ref_seq, in_data.ref,
        SparseVectorInt(static_cast<Eigen::Index>(2 * ref_seq.length())));

    // pairwise alignment for each leaf
    std::string node_seq;
    for(std::size_t node = 0; node < tree.size(); node++) {
        if(tree[node].is_leaf && (tree[node].label != in_data.ref)) {
            float branch = distance_ref(tree, ref_pos, node);
            if(!find_seq(tree[node].label, in_data.fasta_file, node_seq)) {
                throw std::invalid_argument("sequence " + tree[node].label +
                                            " not found in fasta file.");
            }

            pair_seqs[1] = node_seq;

            // P matrix
            if(!in_data.rate.empty()) {
                in_data.mut_model = "user_marg_model";
                P = parse_matrix_csv(in_data.rate);
            } else if(in_data.mut_model.compare("m-ecm") == 0) {
                P = ecm_p(branch, in_data.omega);
            } else {  // m-coati
                P = mg94_p(branch, in_data.omega);
            }

            aln_tmp.f.seq_data.clear();
            if(mg94_marginal(pair_seqs, aln_tmp, P, in_data.frameshifts) != 0) {
                std::cout << "Error: aligning reference " << in_data.ref
                          << " and " << tree[node].label << std::endl;
            }

            SparseVectorInt ins_vector(
                static_cast<Eigen::Index>(aln_tmp.f.seq_data[1].length()));
            insertion_flags(aln_tmp.f.seq_data[0], aln_tmp.f.seq_data[1],
                            ins_vector);

            nodes_ins[node] = insertion_data_t(aln_tmp.f.seq_data[1],
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
    for(const auto& name : in_data.fasta_file.seq_names) {
        auto it = find(nodes_ins[root].names.begin(),
                       nodes_ins[root].names.end(), name);
        auto index = distance(nodes_ins[root].names.begin(), it);
        aln.f.seq_names.push_back(nodes_ins[root].names[index]);
        aln.f.seq_data.push_back(nodes_ins[root].sequences[index]);
    }

    // write alignment
    if(std::filesystem::path(aln.f.path).extension() == ".fasta") {
        return write_fasta(aln.f);
    }
    return write_phylip(aln.f);
}

TEST_CASE("ref_indel_alignment") {
    std::ofstream outfile;
    outfile.open("tree-msa.newick");
    REQUIRE(outfile);
    outfile << "((((A:0.1,B:0.1):0.1,C:0.1):0.1,D:0.1):0.1,E:0.1);";
    outfile.close();

    SUBCASE("m-coati model") {
        input_t input_data("", {"A", "B", "C", "D", "E"},
                           {"TCATCG", "TCAGTCG", "TATCG", "TCACTCG", "TCATC"},
                           "m-coati", "score.log", "test-mcoati-msa.fasta",
                           "tree-msa.newick", "A");
        fasta_t result(input_data.out_file);

        REQUIRE(ref_indel_alignment(input_data) == 0);
        REQUIRE(read_fasta(result) == 0);

        CHECK(std::filesystem::remove(input_data.out_file));

        CHECK(result.seq_names[0] == "A");
        CHECK(result.seq_names[1] == "B");
        CHECK(result.seq_names[2] == "C");
        CHECK(result.seq_names[3] == "D");
        CHECK(result.seq_names[4] == "E");

        CHECK(result.seq_data[0] == "TCA--TCG");
        CHECK(result.seq_data[1] == "TCA-GTCG");
        CHECK(result.seq_data[2] == "T-A--TCG");
        CHECK(result.seq_data[3] == "TCAC-TCG");
        CHECK(result.seq_data[4] == "TCA--TC-");
    }

    SUBCASE("m-ecm model") {
        input_t input_data("", {"A", "B", "C", "D", "E"},
                           {"TCATCG", "TCAGTCG", "TATCG", "TCACTCG", "TCATC"},
                           "m-ecm", "score.log", "test-mecm-msa.fasta",
                           "tree-msa.newick", "A");
        fasta_t result(input_data.out_file);

        REQUIRE(ref_indel_alignment(input_data) == 0);
        REQUIRE(read_fasta(result) == 0);

        CHECK(std::filesystem::remove(input_data.out_file));

        CHECK(result.seq_names[0] == "A");
        CHECK(result.seq_names[1] == "B");
        CHECK(result.seq_names[2] == "C");
        CHECK(result.seq_names[3] == "D");
        CHECK(result.seq_names[4] == "E");

        CHECK(result.seq_data[0] == "TCA--TCG");
        CHECK(result.seq_data[1] == "TCA-GTCG");
        CHECK(result.seq_data[2] == "T-A--TCG");
        CHECK(result.seq_data[3] == "TCAC-TCG");
        CHECK(result.seq_data[4] == "TCA--TC-");
    }
}

float alignment_score(std::vector<std::string> alignment, Matrix& P,
                      float gap_open, float gap_extend) {
    if(alignment[0].length() != alignment[1].length()) {
        throw std::invalid_argument(
            "For alignment scoring both sequences must have equal length.");
    }

    int state = 0;
    float weight = 0.0;
    std::string codon;

    float insertion = gap_open;
    float deletion = gap_open;
    float insertion_ext = gap_extend;
    float deletion_ext = gap_extend;

    // P matrix for marginal Muse and Gaut codon model
    Tensor p = mg94_marginal_p(P);

    std::string seq1 = alignment[0];
    boost::erase_all(seq1, "-");
    int gap_n = 0;

    float nuc_freqs[5] = {0.308, 0.185, 0.199, 0.308, 0.25};

    for(int i = 0, aln_len = static_cast<int>(alignment[0].length());
        i < aln_len; i++) {
        codon = seq1.substr(((i - gap_n) / 3) * 3, 3);  // current codon
        switch(state) {
        case 0:
            if(alignment[0][i] == '-') {
                // insertion;
                unsigned char pos = alignment[1][i];
                weight = weight - logf(insertion) -
                         logf(nuc_freqs[nt4_table[pos]]) -
                         logf(1.0f - insertion_ext);
                state = 1;
                gap_n++;
            } else if(alignment[1][i] == '-') {
                // deletion;
                weight = weight - logf(1.0f - insertion) - logf(deletion) -
                         logf(1.0f - deletion_ext);
                state = 2;
            } else {
                // match/mismatch;
                weight = weight - logf(1.0f - insertion) -
                         logf(1.0f - deletion) -
                         logf(transition(codon, (i + 1 - gap_n) % 3,
                                         alignment[1][i], p));
            }
            break;

        case 1:
            if(alignment[0][i] == '-') {
                // insertion_ext
                unsigned char pos = alignment[1][i];
                weight = weight - logf(insertion_ext) -
                         logf(nuc_freqs[nt4_table[pos]]);
                gap_n++;
            } else if(alignment[1][i] == '-') {
                // deletion
                weight = weight - logf(deletion) - logf(1.0f - deletion_ext);
                state = 2;
            } else {
                // match/mismatch
                weight = weight - logf(1.0f - deletion) -
                         logf(transition(codon, (i + 1 - gap_n) % 3,
                                         alignment[1][i], p));
                state = 0;
            }
            break;

        case 2:
            if(alignment[0][i] == '-') {
                throw std::runtime_error(
                    "Insertion after deletion is not modeled.");
            } else if(alignment[1][i] == '-') {
                // deletion_ext
                weight = weight - logf(deletion_ext);
            } else {
                // match/mismatch
                weight = weight - logf(transition(codon, (i + 1 - gap_n) % 3,
                                                  alignment[1][i], p));
                state = 0;
            }
        }
    }

    return (weight);
}

TEST_CASE("alignment_score") {
    Matrix P(mg94_p(0.0133, 0.2));

    REQUIRE(alignment_score({"CTCTGGATAGTG", "CT----ATAGTG"}, P, 0.001,
                            1.f - 1.f / 6.f) == doctest::Approx(9.29064));
    REQUIRE(alignment_score({"CTCT--AT", "CTCTGGAT"}, P, 0.001,
                            1.f - 1.f / 6.f) == doctest::Approx(12.1493));

    REQUIRE_THROWS_AS(alignment_score({"CTC", "CT"}, P, 0.001, 1.f - 1.f / 6.f),
                      std::invalid_argument);
}
