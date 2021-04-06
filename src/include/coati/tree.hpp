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

#ifndef TREE_H
#define TREE_H

#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <boost/algorithm/string.hpp>

using namespace std;

struct node {
	Eigen::MatrixXd profile;
	vector<string> aln, aln_names;
	bool is_leaf;
	double dist_parent;
	int level;

	node* parent;
	node* lchild;
	node* rchild;

	void set_aln(vector<string> alignment, vector<string> names) {
		aln = alignment;
		aln_names = names;
	}
};

struct tree {
	int leafs, internal_nodes;
	// TODO: think about how to store node names, structure, and edge lengths.
};

int read_newick(string tree_file, tree& guide_tree);

#endif
