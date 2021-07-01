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

using namespace std;
using namespace fst;

/* Alignment using dynamic programming implementation of marginal COATi model */
int mcoati(input_t& in_data, Matrix64f& P) {
    ofstream out_w;
    alignment_t aln;
    aln.f.seq_names = in_data.fasta_file.seq_names;
    aln.f.path = in_data.out_file;

    if(in_data.score) {
        cout << alignment_score(in_data.fasta_file.seq_data, P) << endl;
        return EXIT_SUCCESS;
    }

    if(in_data.mut_model.compare("no_frameshifts") == 0) {
        if(gotoh_noframeshifts(in_data.fasta_file.seq_data, aln, P) != 0) {
            return EXIT_FAILURE;
        }
    } else {
        if(mg94_marginal(in_data.fasta_file.seq_data, aln, P) != 0) {
            return EXIT_FAILURE;
        }
    }

    if(!in_data.weight_file.empty()) {
        // append weight and fasta file name to file
        out_w.open(in_data.weight_file, ios::app | ios::out);
        out_w << in_data.fasta_file.path << "," << in_data.mut_model << ","
              << aln.weight << endl;
        out_w.close();
    }

    // write alignment
    if(aln.f.path.extension() == ".fasta") {
        return write_fasta(aln.f);
    } else {
        return write_phylip(aln.f);
    }
}

TEST_CASE("[align.cc] mcoati") {
    input_t input_data;
    input_data.br_len = 0.0133;
    input_data.fasta_file.seq_names = {"1", "2"};
    input_data.fasta_file.seq_data = {"CTCTGGATAGTG", "CTATAGTG"};
    fasta_t result;
    Matrix64f P;
    mg94_p(P, input_data.br_len);

    SUBCASE("Alignment with frameshifts (default) - output fasta") {
        input_data.out_file = "test-mcoati-fasta.fasta";
        input_data.mut_model = "m-coati";
        input_data.weight_file = "score.log";
        result.path = input_data.out_file;

        if(std::filesystem::exists(input_data.out_file))
            std::filesystem::remove(input_data.out_file);
        if(std::filesystem::exists(input_data.weight_file))
            std::filesystem::remove(input_data.weight_file);

        REQUIRE(mcoati(input_data, P) == 0);
        REQUIRE(read_fasta(result) == 0);

        CHECK(std::filesystem::remove(input_data.out_file));

        CHECK(result.seq_names[0] == "1");
        CHECK(result.seq_names[1] == "2");

        CHECK(result.seq_data[0] == "CTCTGGATAGTG");
        CHECK(result.seq_data[1] == "CT----ATAGTG");

        ifstream infile(input_data.weight_file);
        string s;
        infile >> s;
        CHECK(std::filesystem::remove("score.log"));
        CHECK(s.substr(s.length() - 7) == "9.29064");
    }

    SUBCASE("Alignment with frameshifts (default) - output phylip") {
        input_data.out_file = "test-mcoati-phylip.phy";
        input_data.mut_model = "m-coati";

        if(std::filesystem::exists(input_data.out_file))
            std::filesystem::remove(input_data.out_file);

        REQUIRE(mcoati(input_data, P) == 0);

        ifstream infile(input_data.out_file);
        string s1, s2;

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
    }

    SUBCASE("Alignment with no frameshifts") {
        input_data.fasta_file.seq_data = {"GCGATTGCTGTT", "GCGACTGTT"};
        input_data.out_file = "test-mcoati-no-frameshifts.fasta";
        input_data.mut_model = "no_frameshifts";
        result.path = input_data.out_file;

        if(std::filesystem::exists(input_data.out_file))
            std::filesystem::remove(input_data.out_file);

        REQUIRE(mcoati(input_data, P) == 0);
        REQUIRE(read_fasta(result) == 0);

        CHECK(std::filesystem::remove(input_data.out_file));

        CHECK(result.seq_names[0] == "1");
        CHECK(result.seq_names[1] == "2");

        CHECK(result.seq_data[0] == "GCGATTGCTGTT");
        CHECK(result.seq_data[1] == "GCGA---CTGTT");
    }

    SUBCASE("Score alignment") {
        input_data.out_file = "test-mcoati-score.fasta";
        input_data.mut_model = "m-coati";
        result.path = input_data.out_file;

        if(std::filesystem::exists(input_data.out_file))
            std::filesystem::remove(input_data.out_file);

        REQUIRE(mcoati(input_data, P) == 0);
        REQUIRE(read_fasta(result) == 0);

        CHECK(std::filesystem::remove(input_data.out_file));

        CHECK(alignment_score(result.seq_data, P) == doctest::Approx(9.29064));
    }
}

/* Alignment using FST library*/
int fst_alignment(input_t& in_data, vector<VectorFstStdArc>& fsts) {
    VectorFstStdArc mut_fst;

    if(in_data.mut_model.compare("coati") == 0) {
        mg94(mut_fst, in_data.br_len);
    } else if(in_data.mut_model.compare("dna") == 0) {
        dna(mut_fst, in_data.br_len);
    } else if(in_data.mut_model.compare("ecm") == 0) {
        ecm(mut_fst, in_data.br_len);
    } else {
        cout << "Mutation model unknown. Exiting!" << endl;
        return EXIT_FAILURE;
    }

    // get indel FST
    VectorFstStdArc indel_fst;
    indel(indel_fst, in_data.mut_model);

    // sort mutation and indel FSTs
    VectorFstStdArc mutation_sort, indel_sort;
    mutation_sort = ArcSortFst<StdArc, OLabelCompare<StdArc>>(
        mut_fst, OLabelCompare<StdArc>());
    indel_sort = ArcSortFst<StdArc, ILabelCompare<StdArc>>(
        indel_fst, ILabelCompare<StdArc>());

    // compose mutation and indel FSTs
    ComposeFst<StdArc> coati_comp =
        ComposeFst<StdArc>(mutation_sort, indel_sort);

    // optimize coati FST
    VectorFstStdArc coati_fst;
    coati_fst = optimize(VectorFstStdArc(coati_comp));

    VectorFstStdArc coati_rmep;
    coati_rmep = RmEpsilonFst<StdArc>(coati_fst);  // epsilon removal

    // find alignment graph
    // 1. compose in_tape and coati FSTs
    ComposeFst<StdArc> aln_inter = ComposeFst<StdArc>(fsts[0], coati_rmep);
    // 2. sort intermediate composition
    VectorFstStdArc aln_inter_sort;
    aln_inter_sort = ArcSortFst<StdArc, OLabelCompare<StdArc>>(
        aln_inter, OLabelCompare<StdArc>());
    // 3. compose intermediate and out_tape FSTs
    VectorFstStdArc graph_fst;
    Compose(aln_inter_sort, fsts[1], &graph_fst);

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
    ShortestPath(graph_fst, &aln_path);

    // shortestdistance = weight of shortestpath
    if(!in_data.weight_file.empty()) {
        vector<StdArc::Weight> distance;
        ofstream out_w;

        ShortestDistance(aln_path, &distance);
        // append weight and fasta file name info in file
        out_w.open(in_data.weight_file, ios::app | ios::out);
        out_w << in_data.fasta_file.path << "," << in_data.mut_model << ","
              << distance[0] << endl;
        out_w.close();
    }

    // topsort path FST
    TopSort(&aln_path);

    fasta_t out_fasta(in_data.out_file, in_data.fasta_file.seq_names);

    // write alignment
    if(out_fasta.path.extension() == ".fasta") {
        return write_fasta(aln_path, out_fasta);
    }
    return write_phylip(aln_path, out_fasta);
}

TEST_CASE("[align.cc] fst_alignment") {
    vector<VectorFstStdArc> fsts;
    VectorFstStdArc fsa0, fsa1;
    input_t input_data;
    input_data.br_len = 0.0133;
    input_data.fasta_file.path = "test-fst-alignment.fasta";
    input_data.fasta_file.seq_names = {"1", "2"};
    input_data.fasta_file.seq_data = {"CTCTGGATAGTG", "CTATAGTG"};
    input_data.out_file = "test-fst-alignment.fasta";
    input_data.weight_file = "score.log";

    CHECK(acceptor("CTCTGGATAGTG", fsa0));
    fsts.push_back(fsa0);
    CHECK(acceptor("CTATAGTG", fsa1));
    fsts.push_back(fsa1);

    fasta_t result;
    result.path = input_data.out_file;

    if(std::filesystem::exists(input_data.out_file))
        std::filesystem::remove(input_data.out_file);
    if(std::filesystem::exists(input_data.weight_file))
        std::filesystem::remove(input_data.weight_file);

    SUBCASE("coati model, output fasta") {
        input_data.mut_model = "coati";

        REQUIRE(fst_alignment(input_data, fsts) == 0);
        REQUIRE(read_fasta(result) == 0);

        CHECK(std::filesystem::remove(input_data.out_file));

        CHECK(result.seq_names[0] == "1");
        CHECK(result.seq_names[1] == "2");

        CHECK(result.seq_data[0] == "CTCTGGATAGTG");
        CHECK(result.seq_data[1] == "CT----ATAGTG");

        ifstream infile(input_data.weight_file);
        string s;
        infile >> s;
        CHECK(std::filesystem::remove(input_data.weight_file));
        CHECK(s.substr(s.length() - 7) == "9.31397");
    }

    SUBCASE("coati model, output phylip") {
        input_data.mut_model = "coati";
        input_data.out_file = "test-fst-phylip.phy";

        REQUIRE(fst_alignment(input_data, fsts) == 0);

        ifstream infile(input_data.out_file);
        string s1, s2;

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
        input_data.mut_model = "dna";

        REQUIRE(fst_alignment(input_data, fsts) == 0);
        REQUIRE(read_fasta(result) == 0);

        CHECK(std::filesystem::remove(input_data.out_file));

        CHECK(result.seq_names[0] == "1");
        CHECK(result.seq_names[1] == "2");

        CHECK(result.seq_data[0] == "CTCTGGATAGTG");
        CHECK(result.seq_data[1] == "CT----ATAGTG");

        ifstream infile(input_data.weight_file);
        string s;
        infile >> s;
        CHECK(std::filesystem::remove(input_data.weight_file));
        CHECK(s.substr(s.length() - 7) == "9.31994");
    }

    SUBCASE("ecm model") {
        input_data.mut_model = "ecm";

        REQUIRE(fst_alignment(input_data, fsts) == 0);
        REQUIRE(read_fasta(result) == 0);

        CHECK(std::filesystem::remove(input_data.out_file));

        CHECK(result.seq_names[0] == "1");
        CHECK(result.seq_names[1] == "2");

        CHECK(result.seq_data[0] == "CTCTGGATAGTG");
        CHECK(result.seq_data[1] == "CT----ATAGTG");

        ifstream infile(input_data.weight_file);
        string s;
        infile >> s;
        CHECK(std::filesystem::remove(input_data.weight_file));
        CHECK(s.substr(s.length() - 7) == "9.31388");
    }

    SUBCASE("Unknown model") {
        input_data.mut_model = "unknown";

        REQUIRE(fst_alignment(input_data, fsts) == EXIT_FAILURE);
    }
}

/* Initial msa by collapsing indels after pairwise aln with reference */
int ref_indel_alignment(input_t& in_data) {
    Matrix64f P;
    tree_t tree;
    string newick;
    alignment_t aln, aln_tmp;

    aln.f = fasta_t(in_data.out_file);

    // reack newick tree file
    if(!read_newick(in_data.tree, newick)) {
        cout << "Error: reading newick tree failed." << endl;
        exit(EXIT_FAILURE);
    }

    // parse tree into tree_t (vector<node_t>) variable
    if(parse_newick(newick, tree) != 0) {
        cout << "Error: parsing newick tree failed." << endl;
        exit(EXIT_FAILURE);
    }

    // reroot tree
    if(!reroot(tree, in_data.ref)) {
        cout << "Error: re-rooting tree failed." << endl;
        exit(EXIT_FAILURE);
    }

    // find position of ref in tree
    int ref_pos;
    if(!find_node(tree, in_data.ref, ref_pos)) {
        cout << "Error: reference node not found in tree." << endl;
        exit(EXIT_FAILURE);
    }

    // find sequence of ref in in_data
    vector<string> pair_seqs;
    string ref_seq;
    if(!find_seq(in_data.ref, in_data.fasta_file, ref_seq)) {
        cout << "Error: reference sequence " << in_data.ref
             << " not found in fasta file. " << endl;
        exit(EXIT_FAILURE);
    }

    pair_seqs.push_back(ref_seq);
    pair_seqs.push_back(ref_seq);

    // vector to store insertion_data_t for each node in tree
    vector<insertion_data_t> nodes_ins(tree.size());

    // add insertion_data for REF
    nodes_ins[ref_pos] = insertion_data_t(
        ref_seq, in_data.ref, SparseVectorInt(2 * ref_seq.length()));

    // pairwise alignment for each leaf
    string node_seq;
    for(size_t node = 0; node < tree.size(); node++) {
        if(tree[node].is_leaf && (tree[node].label != in_data.ref)) {
            double branch = distance_ref(tree, ref_pos, node);
            if(!find_seq(tree[node].label, in_data.fasta_file, node_seq)) {
                cout << "Error: sequence " << tree[node].label
                     << " not found in fasta file." << endl;
                exit(EXIT_FAILURE);
            }

            pair_seqs[1] = node_seq;

            // P matrix
            if(in_data.mut_model.compare("m-ecm") == 0) {
                ecm_p(P, branch);
            } else {  // m-coati
                mg94_p(P, branch);
            }

            aln_tmp.f.seq_data.clear();
            if(mg94_marginal(pair_seqs, aln_tmp, P) != 0) {
                cout << "Error: aligning reference " << in_data.ref << " and "
                     << tree[node].label << endl;
            }

            SparseVectorInt ins_vector(2 * aln_tmp.f.seq_data[1].length());
            insertion_flags(aln_tmp.f.seq_data[0], aln_tmp.f.seq_data[1],
                            ins_vector);

            nodes_ins[node] = insertion_data_t(aln_tmp.f.seq_data[1],
                                               tree[node].label, ins_vector);
        }
    }

    // get position of inodes in tree and set leafs as visited (true)
    vector<int> inode_indexes;
    std::vector<int> visited(tree.size(), false);  // list of visited nodes

    for(size_t node = 0; node < tree.size(); node++) {
        if(!tree[node].is_leaf)
            inode_indexes.push_back(node);  // add inode position to vector
        else
            visited[node] = true;  // set leafs to visited
    }

    // fill list of children
    for(size_t i = 0; i < tree.size(); i++) {
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

            if(!children_visited)
                continue;  // if all childen of inode have been visited

            visited[inode_pos] = true;

            // if inode only has a child pass information up
            if(tree[inode_pos].children.size() == 1) {
                nodes_ins[inode_pos] = nodes_ins[tree[inode_pos].children[0]];
                continue;
            }

            // create vector of insertion_data_t with children
            vector<insertion_data_t> tmp_ins_data(
                tree[inode_pos].children.size());
            for(size_t i = 0; i < tree[inode_pos].children.size(); i++) {
                tmp_ins_data[i] = nodes_ins[tree[inode_pos].children[i]];
            }

            // run merge_indels(children_ins_data, nodes_ins[inode_pos]);
            nodes_ins[inode_pos] = insertion_data_t();
            merge_indels(tmp_ins_data, nodes_ins[inode_pos]);
        }
    }

    // transfer result data nodes_ins[ROOT] --> aln && order sequences
    int root = tree[ref_pos].parent;
    for(auto name : in_data.fasta_file.seq_names) {
        vector<string>::iterator it = find(nodes_ins[root].names.begin(),
                                           nodes_ins[root].names.end(), name);
        int index = distance(nodes_ins[root].names.begin(), it);
        aln.f.seq_names.push_back(nodes_ins[root].names[index]);
        aln.f.seq_data.push_back(nodes_ins[root].sequences[index]);
    }

    // write alignment
    if(std::filesystem::path(aln.f.path).extension() == ".fasta") {
        return write_fasta(aln.f);
    } else {
        return write_phylip(aln.f);
    }
}

TEST_CASE("[align.cc] ref_indel_alignment") {
    input_t input_data;
    fasta_t result;
    ofstream outfile;

    input_data.fasta_file.seq_names = {"A", "B", "C", "D", "E"};
    input_data.fasta_file.seq_data = {"TCATCG", "TCAGTCG", "TATCG", "TCACTCG",
                                      "TCATC"};
    input_data.tree = "tree-msa.newick";
    input_data.ref = "A";

    outfile.open(input_data.tree);
    REQUIRE(outfile);
    outfile << "((((A:1,B:1):1,C:1):1,D:1):1,E:1);";
    outfile.close();

    SUBCASE("m-coati model") {
        input_data.mut_model = "m-coati";
        input_data.out_file = "example-mcoati-msa.fasta";
        result.path = input_data.out_file;

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
        input_data.mut_model = "m-ecm";
        input_data.out_file = "example-mecm-msa.fasta";
        result.path = input_data.out_file;

        REQUIRE(ref_indel_alignment(input_data) == 0);
        REQUIRE(read_fasta(result) == 0);

        CHECK(std::filesystem::remove(input_data.out_file));

        CHECK(result.seq_names[0] == "A");
        CHECK(result.seq_names[1] == "B");
        CHECK(result.seq_names[2] == "C");
        CHECK(result.seq_names[3] == "D");
        CHECK(result.seq_names[4] == "E");

        CHECK(result.seq_data[0] == "TC--ATCG");
        CHECK(result.seq_data[1] == "TC-AGTCG");
        CHECK(result.seq_data[2] == "T---ATCG");
        CHECK(result.seq_data[3] == "TCA-CTCG");
        CHECK(result.seq_data[4] == "TC--ATC-");
    }
}

float alignment_score(vector<string> alignment, Matrix64f& P) {
    if(alignment[0].length() != alignment[1].length()) {
        cout << "For alignment scoring both sequences must have equal lenght. "
                "Exiting!"
             << endl;
        exit(EXIT_FAILURE);
    }

    int state = 0;
    double weight = 0.0;
    string codon;

    double insertion = 0.001;
    double deletion = 0.001;
    double insertion_ext = 1.0 - (1.0 / 6.0);
    double deletion_ext = 1.0 - (1.0 / 6.0);

    // P matrix for marginal Muse and Gaut codon model
    Eigen::Tensor<double, 3> p(64, 3, 4);
    mg94_marginal_p(p, P);

    string seq1 = alignment[0];
    boost::erase_all(seq1, "-");
    int gap_n = 0;

    Vector5d nuc_freqs;
    nuc_freqs << 0.308, 0.185, 0.199, 0.308, 0.25;

    for(size_t i = 0; i < alignment[0].length(); i++) {
        codon = seq1.substr(((i - gap_n) / 3) * 3, 3);  // current codon
        switch(state) {
        case 0:
            if(alignment[0][i] == '-') {
                // insertion;
                unsigned char pos = alignment[1][i];
                weight = weight - log(insertion) -
                         log(nuc_freqs[nt4_table[pos]]) -
                         log(1.0 - insertion_ext);
                state = 1;
                gap_n++;
            } else if(alignment[1][i] == '-') {
                // deletion;
                weight = weight - log(1.0 - insertion) - log(deletion) -
                         log(1.0 - deletion_ext);
                state = 2;
            } else {
                // match/mismatch;
                weight = weight - log(1.0 - insertion) - log(1.0 - deletion) -
                         log(transition(codon, (i + 1 - gap_n) % 3,
                                        alignment[1][i], p));
            }
            break;

        case 1:
            if(alignment[0][i] == '-') {
                // insertion_ext
                unsigned char pos = alignment[1][i];
                weight = weight - log(insertion_ext) -
                         log(nuc_freqs[nt4_table[pos]]);
                gap_n++;
            } else if(alignment[1][i] == '-') {
                // deletion
                weight = weight - log(deletion) - log(1.0 - deletion_ext);
                state = 2;
            } else {
                // match/mismatch
                weight = weight - log(1.0 - deletion) -
                         log(transition(codon, (i + 1 - gap_n) % 3,
                                        alignment[1][i], p));
                state = 0;
            }
            break;

        case 2:
            if(alignment[0][i] == '-') {
                cout << "Insertion after deletion is not modeled. Exiting!";
                exit(EXIT_FAILURE);
            } else if(alignment[1][i] == '-') {
                // deletion_ext
                weight = weight - log(deletion_ext);
            } else {
                // match/mismatch
                weight = weight - log(transition(codon, (i + 1 - gap_n) % 3,
                                                 alignment[1][i], p));
                state = 0;
            }
        }
    }

    return (weight);
}
