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
#include <coati/tree.hpp>

namespace newick {
namespace x3 = boost::spirit::x3;
// boost spirit x3 resource:
// https://www.boost.org/doc/libs/1_75_0/libs/spirit/doc/x3/html/index.html
// newick parser resource:
// https://github.com/ultimatesource/mutk/blob/mutk-call/src/lib/newick.cpp

	using x3::float_;
	using x3::lit;
	using x3::_val;
	using x3::_attr;
	using x3::attr;
	using x3::alnum;
	using x3::ascii::char_;
	using boost::fusion::at_c;

	// rule declaration
	x3::rule<class label, string> const label = "label";
	x3::rule<class ilabel, string> const ilabel = "ilabel";
	x3::rule<class length, float> const length = "length";
	x3::rule<class leaf, tree_t> const leaf = "leaf";
	x3::rule<class inode, tree_t> const inode = "inode";
	x3::rule<class node, tree_t> const node = "node";
	x3::rule<class tree, tree_t> const tree = "tree";

	// semantic actions
	auto make_leaf = [](auto& ctx) {
		auto label = at_c<0>(_attr(ctx));		// get label
		auto len = at_c<1>(_attr(ctx));			// get br length
		_val(ctx) = tree_t(1,{label,len,true});	// set to tree_t with 1 leaf node

	};

	auto make_inode = [](auto& ctx) {
		auto label = at_c<1>(_attr(ctx));	// get label
		auto len = at_c<2>(_attr(ctx));		// get br length
		_val(ctx) = tree_t(1,{label,len});	// set to tree_t with 1 (i)node
		auto v = at_c<0>(_attr(ctx));		// get node list (tree_t)
		for(auto &&w : v) {
			auto n = _val(ctx).size();
			for(auto &&x : w) {
				_val(ctx).push_back(x);
				_val(ctx).back().parent += n;
			}
			_val(ctx)[n].parent = 0;
		}

	};

	// rule definition
	auto const tree_def = node >> -lit(';');
	auto const node_def = leaf | inode;
	auto const leaf_def = (label >> length)[make_leaf];
	auto const inode_def = (('(' >> (node % ',') >> ')') >> ilabel >> length)[make_inode];

	auto const label_def = +char_("-0-9A-Za-z/%_.");

	auto const ilabel_def = label | attr("");
	auto const length_def = (':' >> float_) | attr(0.0);

	BOOST_SPIRIT_DEFINE(label);
	BOOST_SPIRIT_DEFINE(ilabel);
	BOOST_SPIRIT_DEFINE(length);
	BOOST_SPIRIT_DEFINE(leaf);
	BOOST_SPIRIT_DEFINE(inode);
	BOOST_SPIRIT_DEFINE(node);
	BOOST_SPIRIT_DEFINE(tree);

} // namespace newick

bool read_newick(string tree_file, string& content) {
	ifstream input(tree_file);		// open input stream
	if(!input.good()) {
		cout << "Error opening '" << tree_file << "'." << endl;
		return false;
	}

	// read newick tree file
	string text((istreambuf_iterator<char>(input)), istreambuf_iterator<char>());
	content = text;

	if(content.length() == 0) {	// Check file isn't empty
		cout << "Reading tree failed, file is empty!" << endl;
		return false;
	}

	return true;
}

/* Read Newick format tree. NOTE: quotation marks on labels not supported.*/
int parse_newick(string content, tree_t& guide_tree) {

	// remove tabs \t ,new lines \n, and spaces
	boost::algorithm::erase_all(content,"\t");
	boost::algorithm::erase_all(content,"\n");
	boost::algorithm::erase_all(content," ");

	auto it = content.begin();
	auto end = content.end();

	bool result = boost::spirit::x3::parse(it, end, newick::tree, guide_tree);

	if(!(result && it == end)) {
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}

TEST_CASE("[tree.cc] parse_newick") {
	tree_t tree;

	REQUIRE(parse_newick("(B_b:6.0,(A-a:5.0,C/c:3.0,E.e:4.0)Ancestor:5.0,D%:11.0);",tree) == 0);
	REQUIRE(tree.size() == 7);

	CHECK(tree[0].label.compare(""));
	CHECK(tree[0].length == 0);
	CHECK(tree[0].is_leaf == false);
	CHECK(tree[0].parent == 0);

	CHECK(tree[1].label.compare("B_b") == 0);
	CHECK(tree[1].length == 6);
	CHECK(tree[1].is_leaf == true);
	CHECK(tree[1].parent == 0);

	CHECK(tree[2].label.compare("Ancestor") == 0);
	CHECK(tree[2].length == 5);
	CHECK(tree[2].is_leaf == false);
	CHECK(tree[2].parent == 0);

	CHECK(tree[3].label.compare("A-a") == 0);
	CHECK(tree[3].length == 5);
	CHECK(tree[3].is_leaf == true);
	CHECK(tree[3].parent == 2);

	CHECK(tree[4].label.compare("C/c") == 0);
	CHECK(tree[4].length == 3);
	CHECK(tree[4].is_leaf == true);
	CHECK(tree[4].parent == 2);

	CHECK(tree[5].label.compare("E.e") == 0);
	CHECK(tree[5].length == 4);
	CHECK(tree[5].is_leaf == true);
	CHECK(tree[5].parent == 2);

	CHECK(tree[6].label.compare("D%") == 0);
	CHECK(tree[6].length == 11);
	CHECK(tree[6].is_leaf == true);
	CHECK(tree[6].parent == 0);
}

/* Determine order of leafs for progressive alignment */
int aln_order(tree_t& tree, vector<pair<int,double>>& order_list) {
	// Part1: find closest pair of leafs

	for(int i=1; i<tree.size(); i++) {  // fill list of children
		tree[tree[i].parent].children.push_back(i);
	}

	pair<int,int> closest_pair;
	float d = FLT_MAX;

	for(int i=0; i<tree.size(); i++) {	// for each node in tree
		if(tree[i].children.empty()) continue; // if no descendants skip
		for(int j=0; j<tree[i].children.size()-1; j++) { // look for closest pair
			for(int k=j+1; k<tree[i].children.size(); k++) {
				// if both nodes aren't children skip
				if(!(tree[tree[i].children[j]].is_leaf && tree[tree[i].children[k]].is_leaf)) {
					continue;
				}
				// if distance between nodes is less than d, update closest pair & d
				if(tree[tree[i].children[j]].length+tree[tree[i].children[k]].length < d) {
					d = tree[tree[i].children[j]].length+tree[tree[i].children[k]].length;
					closest_pair = make_pair(tree[i].children[j], tree[i].children[k]);
				}
			}
		}
	}

	order_list.push_back(make_pair(closest_pair.first,0));
	order_list.push_back(make_pair(closest_pair.second,tree[closest_pair.first].length +
		tree[closest_pair.second].length));

	//Part2: determine order of remaining leafs

	bool visited[tree.size()] = {false}; // list of visited nodes
	visited[order_list[0].first] = visited[order_list[1].first] = true;
	int ancestor = tree[order_list.back().first].parent;
	double branch = 0;

	// while not all nodes have been visited (any value in visitied is false)
	while(any_of(visited, visited+tree.size(), [](bool b){return !b;})) {
		// find leafs
		for(int i=0; i<tree[ancestor].children.size(); i++) {
			// if children is not visited and it's leaf, visit
			if(!visited[tree[ancestor].children[i]] && tree[tree[ancestor].children[i]].is_leaf) {
				visited[tree[ancestor].children[i]] = true;

				order_list.push_back(make_pair(tree[ancestor].children[i],
					tree[tree[ancestor].children[i]].length + branch));
				branch = 0;
			}
		}

		// if any children is inode (not visited on foor lop above) go down that branch
		if(any_of(tree[ancestor].children.begin(),tree[ancestor].children.end(),
			[&visited](int c){return !visited[c];})) {

			for(auto n: tree[ancestor].children) {	// for each children
				if(!visited[n] && !tree[n].is_leaf) { // if !visited & inode, go down branch
					ancestor = n;
					visited[n] = true;
					branch += tree[ancestor].length;
					break;
				}
			}

		} else { // else mark ancestor as visited & move to its parent
			visited[ancestor] = true;
			branch += tree[ancestor].length;
			ancestor = tree[ancestor].parent;
		}
	}

	return EXIT_SUCCESS;
}

TEST_CASE("[tree.cc] aln_order") {
	tree_t tree;
	vector<pair<int,double>> order_list;
	//tree: "(B_b:6.0,(A-a:5.0,C/c:3.0,E.e:4.0)Ancestor:5.0,D%:11.0);"

	tree.push_back(node_t("",0,false));
	tree[0].parent = 0;
	tree.push_back(node_t("B_b",6,true));
	tree[1].parent = 0;
	tree.push_back(node_t("Ancestor",5,false));
	tree[2].parent = 0;
	tree.push_back(node_t("A-a",5,true));
	tree[3].parent = 2;
	tree.push_back(node_t("C/c",3,true));
	tree[4].parent = 2;
	tree.push_back(node_t("E.e",4,true));
	tree[5].parent = 2;
	tree.push_back(node_t("D%",11,true));
	tree[6].parent = 0;

	REQUIRE(aln_order(tree, order_list) == 0);
	REQUIRE(order_list.size() == 5);
	CHECK(order_list[0].first == 4);
	CHECK(order_list[0].second == 0);
	CHECK(order_list[1].first == 5);
	CHECK(order_list[1].second == 7);
	CHECK(order_list[2].first == 3);
	CHECK(order_list[2].second == 5);
	CHECK(order_list[3].first == 1);
	CHECK(order_list[3].second == 11);
	CHECK(order_list[4].first == 6);
	CHECK(order_list[4].second == 11);
}

/* Find fasta sequence given its name */
bool find_seq(string name, fasta_t& f, string& seq) {
	seq.clear();

	for(int i=0; i < f.seq_names.size(); i++) {
		if(f.seq_names[i].compare(name) ==0) {
			seq = f.seq_data[i];
		}
	}

	return !seq.empty();
}

TEST_CASE("[tree.cc] find_seq") {
	string sequence;
	fasta_t fasta = fasta_t("", {"A","B","C"}, {"ACGT","CGTA","GTAC"});

	REQUIRE(!find_seq("Z", fasta, sequence)); // fails, Z is not found -> seq is empty
	REQUIRE(find_seq("A", fasta, sequence));
	CHECK(sequence.compare("ACGT") == 0);
	REQUIRE(find_seq("B", fasta, sequence));
	CHECK(sequence.compare("CGTA") == 0);
	REQUIRE(find_seq("C", fasta, sequence));
	CHECK(sequence.compare("GTAC") == 0);
}
