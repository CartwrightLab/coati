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
#include <coati/tree.hpp>

/* Read Newick format tree. NOTE: quotation marks not supported.*/
int read_newick(string tree_file, tree& guide_tree) {
	ifstream input(tree_file);		// read file
	if(!input.good()) {
		cerr << "Error opening '" << tree_file << "'." << endl;
		return EXIT_FAILURE;
	}

	// Read newick tree file
	string content((istreambuf_iterator<char>(input)), istreambuf_iterator<char>());

	if(content.length() == 0) {	// Check file isn't empty
		cout << "Reading tree failed, file is empty!" << endl;
		return EXIT_FAILURE;
	}
	// Remove tabs \t and new lines \n
	boost::algorithm::erase_all(content,"\t");
	boost::algorithm::erase_all(content,"\n");

	// scan tree looking for operators of interest: parenthesis and commas.
	int left_p, right_p, internal_nodes, edges = 0;
	int leafs = 1;
	vector<int> key_pos;

	for(int i=0; i < content.length(); i++) {
		switch(content.at(i)) {
			case '(':
				key_pos.push_back(i);
				left_p++;
				break;
			case ',':
				key_pos.push_back(i);
				leafs++;
				break;
			case ')':
				key_pos.push_back(i);
				right_p++;
				internal_nodes++;
				break;

		}
	}

	if(left_p != right_p) {
		cout << "Error in newick tree: number of left and right parentheses does not match"
			<< endl;
		return EXIT_FAILURE;
	}

	edges = leafs + internal_nodes - 1;

	// parsing of Newick tree
	int pos;
	for(int i=0; i < key_pos.size(); i++) {
		pos = key_pos[i];
		cout << "pos: " << pos << endl;
		// TODO: figure out different case scenarios
	}

	return EXIT_SUCCESS;
}
