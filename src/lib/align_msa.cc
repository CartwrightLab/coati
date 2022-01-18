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

#include <coati/align_msa.hpp>

namespace coati {

/**
 * \brief Pairwise alignment using dynamic programming and a marginal model.
 *
 * Multiple sequence alignment by pairwise alignment reference with rest of
 *  sequences using coati::align_pair, then merging indels.
 *
 * @param[in] input coati::alignment_t alignment parameters.
 */
/* Initial msa by collapsing indels after pairwise aln with reference */
bool ref_indel_alignment(coati::alignment_t& input) {
    coati::tree::tree_t tree;
    std::string newick;
    coati::alignment_t aln, aln_tmp;

    aln.data.out_file = input.data.out_file;

    // read newick tree file
    if(!coati::tree::read_newick(input.tree, newick)) {
        throw std::invalid_argument("Reading newick tree failed.");
    }

    // parse tree into tree_t (vector<node_t>) variable
    if(coati::tree::parse_newick(newick, tree) != 0) {
        throw std::invalid_argument("Parsing newick tree failed.");
    }

    // reroot tree
    if(!coati::tree::reroot(tree, input.ref)) {
        throw std::invalid_argument("Re-rooting tree failed.");
    }

    // find position of ref in tree
    std::size_t ref_pos = 0;
    if(!coati::tree::find_node(tree, input.ref, ref_pos)) {
        throw std::invalid_argument("Reference node not found in tree.");
    }

    // find sequence of ref in input
    std::vector<std::string> pair_seqs;
    std::string ref_seq;
    if(!coati::tree::find_seq(input.ref, input.data, ref_seq)) {
        throw std::invalid_argument("reference sequence " + input.ref +
                                    " not found in input file.");
    }

    pair_seqs.push_back(ref_seq);
    pair_seqs.push_back(ref_seq);

    // vector to store insertion_data_t for each node in tree
    std::vector<insertion_data_t> nodes_ins(tree.size());

    // add insertion_data for REF
    nodes_ins[ref_pos] = insertion_data_t(
        ref_seq, input.ref,
        SparseVectorInt(static_cast<Eigen::Index>(2 * ref_seq.length())));

    // pairwise alignment for each leaf
    std::string node_seq;
    for(std::size_t node = 0; node < tree.size(); node++) {
        if(tree[node].is_leaf && (tree[node].label != input.ref)) {
            float branch = distance_ref(tree, ref_pos, node);
            if(!coati::tree::find_seq(tree[node].label, input.data, node_seq)) {
                throw std::invalid_argument("sequence " + tree[node].label +
                                            " not found in input file.");
            }

            pair_seqs[1] = node_seq;

            input.br_len = branch;
            coati::utils::set_subst(input);

            aln_tmp.data.seqs.clear();
            auto anc = pair_seqs[0];
            auto des = pair_seqs[1];
            coati::utils::sequence_pair_t seq_pair =
                coati::utils::marginal_seq_encoding(anc, des);
            coati::align_pair_work_mem_t work;
            coati::viterbi_mem(work, seq_pair[0], seq_pair[1], input);
            coati::traceback(work, anc, des, aln_tmp, input.gap.len);

            SparseVectorInt ins_vector(
                static_cast<Eigen::Index>(aln_tmp.data.seqs[1].length()));
            insertion_flags(aln_tmp.data.seqs[0], aln_tmp.data.seqs[1],
                            ins_vector);

            nodes_ins[node] = insertion_data_t(aln_tmp.data.seqs[1],
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
    for(const auto& name : input.data.names) {
        auto it = find(nodes_ins[root].names.begin(),
                       nodes_ins[root].names.end(), name);
        auto index = distance(nodes_ins[root].names.begin(), it);
        aln.data.names.push_back(nodes_ins[root].names[index]);
        aln.data.seqs.push_back(nodes_ins[root].sequences[index]);
    }

    // write alignment
    return coati::utils::write_output(aln.data);
}

/// @private
TEST_CASE("ref_indel_alignment") {
    std::ofstream outfile;
    outfile.open("tree-msa.newick");
    REQUIRE(outfile);
    outfile << "((((A:0.1,B:0.1):0.1,C:0.1):0.1,D:0.1):0.1,E:0.1);";
    outfile.close();

    SUBCASE("m-coati model") {
        coati::alignment_t aln;
        aln.data =
            coati::data_t("", {"A", "B", "C", "D", "E"},
                          {"TCATCG", "TCAGTCG", "TATCG", "TCACTCG", "TCATC"});
        aln.model = "m-coati";
        aln.data.out_file = {{"test-mecm-msa.fasta"}, {".fasta"}};
        aln.tree = "tree-msa.newick";
        aln.ref = "A";

        REQUIRE(ref_indel_alignment(aln));

        std::ifstream infile(aln.data.out_file.path);
        std::string s1, s2;

        infile >> s1 >> s2;
        CHECK(s1 == ">A");
        CHECK(s2 == "TCA--TCG");
        infile >> s1 >> s2;
        CHECK(s1 == ">B");
        CHECK(s2 == "TCA-GTCG");
        infile >> s1 >> s2;
        CHECK(s1 == ">C");
        CHECK(s2 == "T-A--TCG");
        infile >> s1 >> s2;
        CHECK(s1 == ">D");
        CHECK(s2 == "TCAC-TCG");
        infile >> s1 >> s2;
        CHECK(s1 == ">E");
        CHECK(s2 == "TCA--TC-");
        CHECK(std::filesystem::remove(aln.data.out_file.path));
    }

    SUBCASE("m-ecm model") {
        coati::alignment_t aln;
        aln.data =
            coati::data_t("", {"A", "B", "C", "D", "E"},
                          {"TCATCG", "TCAGTCG", "TATCG", "TCACTCG", "TCATC"});
        aln.model = "m-ecm";
        aln.data.out_file = {{"test-mecm-msa.fasta"}, {".fasta"}};
        aln.tree = "tree-msa.newick";
        aln.ref = "A";

        REQUIRE(ref_indel_alignment(aln));

        std::ifstream infile(aln.data.out_file.path);
        std::string s1, s2;

        infile >> s1 >> s2;
        CHECK(s1 == ">A");
        CHECK(s2 == "TCA--TCG");
        infile >> s1 >> s2;
        CHECK(s1 == ">B");
        CHECK(s2 == "TCA-GTCG");
        infile >> s1 >> s2;
        CHECK(s1 == ">C");
        CHECK(s2 == "T-A--TCG");
        infile >> s1 >> s2;
        CHECK(s1 == ">D");
        CHECK(s2 == "TCAC-TCG");
        infile >> s1 >> s2;
        CHECK(s1 == ">E");
        CHECK(s2 == "TCA--TC-");
        CHECK(std::filesystem::remove(aln.data.out_file.path));
    }

    CHECK(std::filesystem::remove("tree-msa.newick"));
}
}  // namespace coati
