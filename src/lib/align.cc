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
    if(boost::filesystem::extension(aln.f.path) == ".fasta") {
        return write_fasta(aln.f);
    } else {
        return write_phylip(aln.f);
    }
}

/* Alignment using FST library*/
int fst_alignment(input_t& in_data, vector<VectorFst<StdArc>>& fsts) {
    VectorFst<StdArc> mut_fst;

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
    VectorFst<StdArc> indel_fst;
    indel(indel_fst, in_data.mut_model);

    // sort mutation and indel FSTs
    VectorFst<StdArc> mutation_sort, indel_sort;
    mutation_sort = ArcSortFst<StdArc, OLabelCompare<StdArc>>(
        mut_fst, OLabelCompare<StdArc>());
    indel_sort = ArcSortFst<StdArc, ILabelCompare<StdArc>>(
        indel_fst, ILabelCompare<StdArc>());

    // compose mutation and indel FSTs
    ComposeFst<StdArc> coati_comp =
        ComposeFst<StdArc>(mutation_sort, indel_sort);

    // optimize coati FST
    VectorFst<StdArc> coati_fst;
    coati_fst = optimize(VectorFst<StdArc>(coati_comp));

    // find alignment graph
    // 1. compose in_tape and coati FSTs
    ComposeFst<StdArc> aln_inter = ComposeFst<StdArc>(fsts[0], coati_fst);
    // 2. sort intermediate composition
    VectorFst<StdArc> aln_inter_sort;
    aln_inter_sort = ArcSortFst<StdArc, OLabelCompare<StdArc>>(
        aln_inter, OLabelCompare<StdArc>());
    // 3. compose intermediate and out_tape FSTs
    VectorFst<StdArc> graph_fst;
    Compose(aln_inter_sort, fsts[1], &graph_fst);

    //  case 1: ComposeFst is time: O(v1 v2 d1 (log d2 + m2)), space O(v1 v2)
    //			ShortestPath with ComposeFst (PDT) is time: O((V+E)^3),
    //space: O((V+E)^3)
    //          Then convert to VectorFst (FST) is time, space: O(e^(O(V+E)))
    //  case 2: ComposeFst is time: O(v1 v2 d1 (log d2 + m2)), space: O(v1 v2)
    //			Convert ComposeFst (PDT) to VectorFst(FST) is time,
    //space:
    // O(e^(O(V+E))) 			Then ShortestPath with VectorFst is O(V log(V) +
    // E)
    //  case 3: Compose is time: O(V1 V2 D1 (log D2 + M2)), space O(V1 V2 D1 M2)
    //			Then ShortestPath with VectorFst is O(V log(V) + E)

    // find shortest path through graph
    VectorFst<StdArc> aln_path;
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
    if(boost::filesystem::extension(out_fasta.path) == ".fasta") {
        return write_fasta(aln_path, out_fasta);
    } else {
        return write_phylip(aln_path, out_fasta);
    }
}

/* Progressive alignment  */
int progressive_aln(input_t& in_data) {
    Matrix64f P;
    vector<pair<int, double>> order;
    tree_t tree;
    string newick;
    ofstream out_w;
    alignment_t aln;

    aln.f.seq_names = in_data.fasta_file.seq_names;
    aln.f.path = in_data.out_file;

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

    // get order of leafs to align and branch lengths
    aln_order(tree, order);

    // retrieve closest leafs' sequences
    vector<string> closest_leafs, next_leaf;
    string seq;

    // sequence for 1st leaf
    if(!find_seq(tree[order[0].first].label, in_data.fasta_file, seq)) {
        cout << "Error: sequence " << tree[order[0].first].label
             << " not find in fasta file." << endl;
        exit(EXIT_FAILURE);
    }
    closest_leafs.push_back(seq);
    aln.f.seq_names[0] = tree[order[0].first].label;  // label (name) of leaf

    // sequence for 2nd leaf
    if(!find_seq(tree[order[1].first].label, in_data.fasta_file, seq)) {
        cout << "Error: sequence " << tree[order[1].first].label
             << " not find in fasta file." << endl;
        exit(EXIT_FAILURE);
    }
    closest_leafs.push_back(seq);
    aln.f.seq_names[1] = tree[order[1].first].label;  // label (name) of leaf

    // get branch length and create P matrix
    double branch = order[1].second;
    if(in_data.mut_model.compare("m-ecm") == 0) {
        ecm_p(P, branch);
    } else {  // m-coati
        mg94_p(P, branch);
    }

    // alignpair closest leafs
    if(mg94_marginal(closest_leafs, aln, P) != 0) {
        cout << "Error: aligning closest leafs " << tree[order[0].first].label
             << " and " << tree[order[1].first].label << endl;
        exit(EXIT_FAILURE);
    }

    // align rest of leafs in order
    for(int i = 2; i < order.size(); i++) {
        // find sequence of next leaf
        if(!find_seq(tree[order[i].first].label, in_data.fasta_file, seq)) {
            cout << "Error: sequence " << tree[order[i].first].label
                 << " not find in fasta file." << endl;
            exit(EXIT_FAILURE);
        }
        next_leaf = {seq};
        // new branch length and update P matrix
        if(in_data.mut_model.compare("m-ecm") == 0) {
            ecm_p(P, order[i].second);
        } else {  // m-coati
            mg94_p(P, order[i].second);
        }

        // align
        gotoh_profile_marginal(aln.f.seq_data, next_leaf, aln, P);
        aln.f.seq_names[i] =
            tree[order[i].first].label;  // add sequence label (name)
    }

    // if file to save weight of alignment is provided do so
    if(!in_data.weight_file.empty()) {
        // append weight and fasta file name to file
        out_w.open(in_data.weight_file, ios::app | ios::out);
        out_w << in_data.fasta_file.path << "," << in_data.mut_model << ","
              << aln.weight << endl;
        out_w.close();
    }

    // write alignment
    if(boost::filesystem::extension(aln.f.path) == ".fasta") {
        return write_fasta(aln.f);
    } else {
        return write_phylip(aln.f);
    }
}
