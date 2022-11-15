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

#include <coati/align_msa.hpp>

namespace coati {

/**
 * @brief Multiple sequence alignment using an iterative algorithm.
 *
 * Initial multiple sequence alignment (to be refined - in development)
 * by collapsing indels along the tree after pairwise alignment with reference
 * sequence.
 *
 * @details Pairwise align all sequences with the refence. Then, starting with
 * the closest pair of leafs the aligned sequences are merged (all in frame
 * with the reference). Deletions are kept, insertions are merged when part
 * of the same branch, kept sepparate otherwise even when contiguous.
 *
 * @param[in] input coati::alignment_t sequences to align and parameters.
 *
 * @retval true successful run.
 */
bool ref_indel_alignment(coati::alignment_t& input) {
    coati::data_t aligned;

    if(!input.is_marginal()) {
        throw std::invalid_argument("MSA only supports marginal models.");
    }

    // read input data
    input.data = coati::io::read_input(input);
    if(input.data.size() < 3) {
        throw std::invalid_argument("At least three sequences required.");
    }

    aligned.out_file = input.data.out_file;

    // read newick tree file
    std::string newick = coati::tree::read_newick(input.tree);

    // parse tree into tree_t (vector<node_t>) variable
    coati::tree::tree_t tree = coati::tree::parse_newick(newick);

    // reroot tree
    coati::tree::reroot(tree, input.refs);

    // find position of ref in tree
    std::size_t ref_pos = coati::tree::find_node(tree, input.refs);

    // find sequence of ref in input
    std::string ref_seq = coati::tree::find_seq(input.refs, input.data);

    // vector to store insertion_data_t for each node in tree
    coati::insertion_vector nodes_ins(tree.size());

    // initialize insertion_data for REF
    nodes_ins[ref_pos] = insertion_data_t(
        ref_seq, input.refs,
        SparseVectorInt(static_cast<Eigen::Index>(2 * ref_seq.length())));

    // pairwise alignment of each leaf with reference and store insertion info.
    align_leafs(input, tree, ref_pos, ref_seq, nodes_ins);

    // get position of inodes in tree and set leafs as visited (true)
    std::vector<std::size_t> inode_indexes;
    std::vector<bool> visited(tree.size(), false);  // list of visited nodes

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

    // merge pairwise alignments up the tree until root - final MSA
    merge_alignments(visited, tree, nodes_ins, inode_indexes);

    // transfer result data nodes_ins[ROOT] --> aln && order sequences
    auto root = tree[ref_pos].parent;
    for(const auto& name : input.data.names) {
        auto it = find(nodes_ins[root].names.begin(),
                       nodes_ins[root].names.end(), name);
        auto index = distance(nodes_ins[root].names.begin(), it);
        aligned.names.push_back(nodes_ins[root].names[index]);
        aligned.seqs.push_back(nodes_ins[root].sequences[index]);
    }

    // write alignment
    coati::io::write_output(aligned);
    return true;
}

/// @private
// GCOVR_EXCL_START
TEST_CASE("ref_indel_alignment") {
    std::ofstream outfile;
    outfile.open("tree-msa.newick");
    REQUIRE(outfile);
    outfile << "((((A:0.1,B:0.1):0.1,C:0.1):0.1,D:0.1):0.1,E:0.1);";
    outfile.close();

    outfile.open("test-msa.fasta");
    REQUIRE(outfile);
    outfile << ">A\nTCATCG\n>B\nTCAGTCG\n>C\nTATCG\n>D\nTCACTCG\n>"
               "E\nTCATC\n";
    outfile.close();

    SUBCASE("marginal model") {
        coati::alignment_t aln;
        aln.data.path = "test-msa.fasta";
        aln.model = "marginal";
        aln.output = "test-mecm-msa.fasta";
        aln.tree = "tree-msa.newick";
        aln.refs = "A";

        REQUIRE(ref_indel_alignment(aln));

        std::ifstream infile(aln.data.out_file.path);
        std::string s1, s2;

        infile >> s1 >> s2;
        CHECK_EQ(s1, ">A");
        CHECK_EQ(s2, "TCA--TCG");
        infile >> s1 >> s2;
        CHECK_EQ(s1, ">B");
        CHECK_EQ(s2, "TCA-GTCG");
        infile >> s1 >> s2;
        CHECK_EQ(s1, ">C");
        CHECK_EQ(s2, "T-A--TCG");
        infile >> s1 >> s2;
        CHECK_EQ(s1, ">D");
        CHECK_EQ(s2, "TCAC-TCG");
        infile >> s1 >> s2;
        CHECK_EQ(s1, ">E");
        CHECK_EQ(s2, "TCA--TC-");
        REQUIRE(std::filesystem::remove(aln.output));
    }

    SUBCASE("m-ecm model") {
        coati::alignment_t aln;
        aln.data.path = "test-msa.fasta";
        aln.model = "m-ecm";
        aln.output = "test-mecm-msa.fasta";
        aln.tree = "tree-msa.newick";
        aln.refs = "A";

        REQUIRE(ref_indel_alignment(aln));

        std::ifstream infile(aln.data.out_file.path);
        std::string s1, s2;

        infile >> s1 >> s2;
        CHECK_EQ(s1, ">A");
        CHECK_EQ(s2, "TCA--TCG");
        infile >> s1 >> s2;
        CHECK_EQ(s1, ">B");
        CHECK_EQ(s2, "TCA-GTCG");
        infile >> s1 >> s2;
        CHECK_EQ(s1, ">C");
        CHECK_EQ(s2, "T-A--TCG");
        infile >> s1 >> s2;
        CHECK_EQ(s1, ">D");
        CHECK_EQ(s2, "TCAC-TCG");
        infile >> s1 >> s2;
        CHECK_EQ(s1, ">E");
        CHECK_EQ(s2, "TCA--TC-");
        REQUIRE(std::filesystem::remove(aln.output));
    }

    SUBCASE("Less than 3 sequences - fail") {
        std::ofstream outfile;
        outfile.open("tree-msa.newick");
        REQUIRE(outfile);
        outfile << "(A:0.1,B:0.1);\n";
        outfile.close();

        outfile.open("test-msa.fasta");
        REQUIRE(outfile);
        outfile << ">A\nTCATCG\n>B\nTCAGTCG\n";
        outfile.close();

        coati::alignment_t aln;
        aln.data.path = "test-msa.fasta";
        aln.tree = "tree-msa.newick";
        aln.refs = "A";

        CHECK_THROWS_AS(ref_indel_alignment(aln), std::invalid_argument);
    }

    SUBCASE("More complex tree") {
        std::ofstream outfile;
        coati::alignment_t aln;
        outfile.open("tree-msa.newick");
        REQUIRE(outfile);
        outfile << "((A:0.1,B:0.1):0.1,(C:0.1,(D:0.1,E:0.1):0.1):0.1,F:0.1);\n";
        outfile.close();

        outfile.open("test-msa.fasta");
        REQUIRE(outfile);
        outfile << ">A\nTCATCG\n>B\nTCAGTCG\n>C\nTATCG\n>D\nTCACTCG\n>"
                   "E\nTCATC\n>F\nTCATCG";
        outfile.close();

        aln.data.path = "test-msa.fasta";
        aln.tree = "tree-msa.newick";
        aln.refs = "A";
        aln.output = "test-fst-complex-tree.fa";
        REQUIRE(ref_indel_alignment(aln));

        std::ifstream infile(aln.data.out_file.path);
        std::string s1, s2;

        infile >> s1 >> s2;
        CHECK_EQ(s1, ">A");
        CHECK_EQ(s2, "TCA--TCG");
        infile >> s1 >> s2;
        CHECK_EQ(s1, ">B");
        CHECK_EQ(s2, "TCA-GTCG");
        infile >> s1 >> s2;
        CHECK_EQ(s1, ">C");
        CHECK_EQ(s2, "T-A--TCG");
        infile >> s1 >> s2;
        CHECK_EQ(s1, ">D");
        CHECK_EQ(s2, "TCAC-TCG");
        infile >> s1 >> s2;
        CHECK_EQ(s1, ">E");
        CHECK_EQ(s2, "TCA--TC-");
        infile >> s1 >> s2;
        CHECK_EQ(s1, ">F");
        CHECK_EQ(s2, "TCA--TCG");
        REQUIRE(std::filesystem::remove(aln.output));
    }
    SUBCASE("not marginal model - fail") {
        coati::alignment_t aln;
        aln.model = "fst";

        CHECK_THROWS_AS(ref_indel_alignment(aln), std::invalid_argument);
    }

    REQUIRE(std::filesystem::remove("tree-msa.newick"));
    REQUIRE(std::filesystem::remove("test-msa.fasta"));
}
// GCOVR_EXCL_STOP

/**
 * @brief Pairwise alignments of leafs with reference sequence.
 *
 * @details Pairwise align each leaf with the refence and store insertions
 * (w.r.t. reference).
 *
 * @param[in,out] input coati::alignment_t sequences to align and parameters.
 * @param[in] tree coati::tree::tree_t phylogenetic tree.
 * @param[in] ref_pos std::size_t position of reference seq in tree.
 * @param[in] ref_seq std::string reference sequence.
 * @param[in,out] nodes_ins coati::insertion_vector position of insertions for
 * each sequence.
 */
void align_leafs(coati::alignment_t& input, const coati::tree::tree_t& tree,
                 std::size_t ref_pos, const std::string& ref_seq,
                 coati::insertion_vector& nodes_ins) {
    coati::alignment_t aln;
    std::vector<std::string> pair_seqs(2);

    pair_seqs[0] = ref_seq;

    for(std::size_t node = 0; node < tree.size(); ++node) {
        // if node is a leaf and not the reference pairwise align w/ reference
        if(tree[node].is_leaf && (tree[node].label != input.refs)) {
            // find branch length in tree (given by user) and sequence of node
            input.br_len = distance_ref(tree, ref_pos, node);
            pair_seqs[1] = coati::tree::find_seq(tree[node].label, input.data);

            coati::utils::set_subst(input);

            // pairwise alignment of reference and leaf
            aln.data.seqs.clear();
            coati::utils::sequence_pair_t seq_pair =
                coati::utils::marginal_seq_encoding(pair_seqs[0], pair_seqs[1]);
            coati::align_pair_work_mem_t work;
            coati::viterbi_mem(work, seq_pair[0], seq_pair[1], input);
            coati::traceback(work, pair_seqs[0], pair_seqs[1], aln,
                             input.gap.len);

            // store insertion positions and type (open by default)
            SparseVectorInt ins = insertion_flags(aln.seq(0), aln.seq(1));

            nodes_ins[node] =
                insertion_data_t(aln.data.seqs[1], tree[node].label, ins);
        }
    }
}

/**
 * @brief Merge alignments starting from leafs until root.
 *
 * @details Merge insertion positions starting with closest nodes and going
 * up the tree until arriving to the root. Insertions that are in the same
 * position and the character inserted is the same are grouped. When, along the
 * tree, an insertion(s) does not match the inserted nucleotide, these are
 * considred closed and no more insertion can share the same position.
 * This ensures events in different branches are not grouped together.
 *
 * @param[in,out] visited std::vetor<bool> nodes visited.
 * @param[in] tree coati::tree::tree_t phylogenetic tree.
 * @param[in,out] nodes_ins coati::insertion_vector position and type of
 * insertions.
 * @param[in] inode_indexes std::vector<std::size_t> index of internal nodes.
 */
void merge_alignments(std::vector<bool>& visited,
                      const coati::tree::tree_t& tree,
                      coati::insertion_vector& nodes_ins,
                      const std::vector<std::size_t>& inode_indexes) {
    // while not all nodes have been visited (any value in visited is false)
    while(any_of(visited.begin(), visited.end(), [](bool b) { return !b; })) {
        for(auto inode_pos : inode_indexes) {  // for all internal nodes
            if(visited[inode_pos]) {
                continue;
            }
            // if any children is not visited, skip & come back when all are
            if(std::any_of(tree[inode_pos].children.begin(),
                           tree[inode_pos].children.end(),
                           [visited](int c) { return !visited[c]; })) {
                continue;
            }

            visited[inode_pos] = true;

            // if inode only has a child pass information up
            if(tree[inode_pos].children.size() == 1) {
                nodes_ins[inode_pos] = nodes_ins[tree[inode_pos].children[0]];
                continue;
            }

            // create vector of insertion_data_t with children
            coati::insertion_vector tmp_ins_data;
            tmp_ins_data.reserve(tree[inode_pos].children.size());
            for(auto child : tree[inode_pos].children) {
                tmp_ins_data.emplace_back(nodes_ins[child]);
            }

            // merge insertions from all children of inode
            //  and store in self (inode) position
            nodes_ins[inode_pos] = insertion_data_t();
            merge_indels(tmp_ins_data, nodes_ins[inode_pos]);
        }
    }
}

}  // namespace coati
