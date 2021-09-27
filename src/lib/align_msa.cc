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

    CHECK(std::filesystem::remove("tree-msa.newick"));
}
}  // namespace coati
