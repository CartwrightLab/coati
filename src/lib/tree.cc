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

#include <cfloat>
#include <coati/tree.hpp>

namespace newick {
namespace x3 = boost::spirit::x3;
// boost spirit x3 resource:
// https://www.boost.org/doc/libs/1_75_0/libs/spirit/doc/x3/html/index.html
// newick parser resource:
// https://github.com/ultimatesource/mutk/blob/mutk-call/src/lib/newick.cpp

using boost::fusion::at_c;
using x3::_attr;
using x3::_val;
using x3::attr;
using x3::float_;
using x3::lit;
using x3::ascii::char_;

// rule declaration
x3::rule<class label, std::string> const label =  // NOLINT(cert-err58-cpp)
    "label";
x3::rule<class ilabel, std::string> const ilabel =  // NOLINT(cert-err58-cpp)
    "ilabel";
x3::rule<class length, float> const length =  // NOLINT(cert-err58-cpp)
    "length";
x3::rule<class leaf, tree_t> const leaf = "leaf";     // NOLINT(cert-err58-cpp)
x3::rule<class inode, tree_t> const inode = "inode";  // NOLINT(cert-err58-cpp)
x3::rule<class node, tree_t> const node = "node";     // NOLINT(cert-err58-cpp)
x3::rule<class tree, tree_t> const tree = "tree";     // NOLINT(cert-err58-cpp)

// semantic actions
auto const make_leaf = [](auto& ctx) {
    auto label = at_c<0>(_attr(ctx));  // get label
    auto len = at_c<1>(_attr(ctx));    // get br length
    _val(ctx) =
        tree_t(1, {label, len, true});  // set to tree_t with 1 leaf node
};

auto const make_inode = [](auto& ctx) {
    auto label = at_c<1>(_attr(ctx));     // get label
    auto len = at_c<2>(_attr(ctx));       // get br length
    _val(ctx) = tree_t(1, {label, len});  // set to tree_t with 1 (i)node
    auto v = at_c<0>(_attr(ctx));         // get node list (tree_t)
    for(auto&& w : v) {
        auto n = _val(ctx).size();
        for(auto&& x : w) {
            _val(ctx).push_back(x);
            _val(ctx).back().parent += n;
        }
        _val(ctx)[n].parent = 0;
    }
};

// rule definition
auto const tree_def = node >> -lit(';');             // NOLINT(cert-err58-cpp)
auto const node_def = leaf | inode;                  // NOLINT(cert-err58-cpp)
auto const leaf_def = (label >> length)[make_leaf];  // NOLINT(cert-err58-cpp)
auto const inode_def =                               // NOLINT(cert-err58-cpp)
    (('(' >> (node % ',') >> ')') >> ilabel >>
     length)[make_inode];  // NOLINT(cert-err58-cpp)

auto const label_def = +char_("-0-9A-Za-z/%_.");  // NOLINT(cert-err58-cpp)

auto const ilabel_def = label | attr("");             // NOLINT(cert-err58-cpp)
auto const length_def = (':' >> float_) | attr(0.0);  // NOLINT(cert-err58-cpp)

BOOST_SPIRIT_DEFINE(label);   // NOLINT(performance-unnecessary-value-param)
BOOST_SPIRIT_DEFINE(ilabel);  // NOLINT(performance-unnecessary-value-param)
BOOST_SPIRIT_DEFINE(length);  // NOLINT(performance-unnecessary-value-param)
BOOST_SPIRIT_DEFINE(leaf);    // NOLINT(performance-unnecessary-value-param)
// NOLINTNEXTLINE(misc-no-recursion, performance-unnecessary-value-param)
BOOST_SPIRIT_DEFINE(inode);
// NOLINTNEXTLINE(misc-no-recursion, performance-unnecessary-value-param)
BOOST_SPIRIT_DEFINE(node);
BOOST_SPIRIT_DEFINE(tree);  // NOLINT(performance-unnecessary-value-param)

}  // namespace newick

bool read_newick(const std::string& tree_file, std::string& content) {
    std::ifstream input(tree_file);  // open input stream
    if(!input.good()) {
        throw std::invalid_argument("Error opening " + tree_file + ".");
    }

    // read newick tree file
    std::string text((std::istreambuf_iterator<char>(input)),
                     std::istreambuf_iterator<char>());
    content = text;

    if(content.length() == 0) {  // Check file isn't empty
        throw std::invalid_argument("Reading tree failed, file is empty!");
    }

    return true;
}

/* Read Newick format tree. NOTE: quotation marks on labels not supported.*/
int parse_newick(std::string content, tree_t& guide_tree) {
    // remove tabs \t ,new lines \n, and spaces
    boost::algorithm::erase_all(content, "\t");
    boost::algorithm::erase_all(content, "\n");
    boost::algorithm::erase_all(content, " ");

    auto it = content.begin();
    auto end = content.end();

    // NOLINTNEXTLINE(misc-no-recursion)
    bool result = boost::spirit::x3::parse(it, end, newick::tree, guide_tree);

    if(!(result && it == end)) {
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

TEST_CASE("[tree.cc] parse_newick") {
    tree_t tree;

    REQUIRE(
        parse_newick("(B_b:6.0,(A-a:5.0,C/c:3.0,E.e:4.0)Ancestor:5.0,D%:11.0);",
                     tree) == 0);
    REQUIRE(tree.size() == 7);

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
int aln_order(tree_t& tree, std::vector<std::pair<int, float>>& order_list) {
    // Part1: find closest pair of leafs

    for(std::size_t i = 1; i < tree.size(); i++) {  // fill list of children
        tree[tree[i].parent].children.push_back(i);
    }

    std::pair<int, int> closest_pair;
    float d = FLT_MAX;

    for(std::size_t i = 0; i < tree.size(); i++) {  // for each node in tree
        if(tree[i].children.empty()) continue;      // if no descendants skip
        for(std::size_t j = 0; j < tree[i].children.size() - 1;
            j++) {  // look for closest pair
            for(std::size_t k = j + 1; k < tree[i].children.size(); k++) {
                // if both nodes aren't children skip
                if(!(tree[tree[i].children[j]].is_leaf &&
                     tree[tree[i].children[k]].is_leaf)) {
                    continue;
                }
                // if distance between nodes is less than d, update closest pair
                // & d
                if(tree[tree[i].children[j]].length +
                       tree[tree[i].children[k]].length <
                   d) {
                    d = tree[tree[i].children[j]].length +
                        tree[tree[i].children[k]].length;
                    closest_pair = std::make_pair(tree[i].children[j],
                                                  tree[i].children[k]);
                }
            }
        }
    }

    order_list.emplace_back(closest_pair.first, 0);
    order_list.emplace_back(
        closest_pair.second,
        tree[closest_pair.first].length + tree[closest_pair.second].length);

    // Part2: determine order of remaining leafs

    std::vector<int> visited(tree.size(), false);  // list of visited nodes

    visited[order_list[0].first] = visited[order_list[1].first] = true;
    std::size_t ancestor = tree[order_list.back().first].parent;
    float branch = 0;

    // while not all nodes have been visited (any value in visitied is false)
    while(any_of(visited.begin(), visited.end(), [](bool b) { return !b; })) {
        // find leafs
        for(std::size_t i = 0; i < tree[ancestor].children.size(); i++) {
            // if children is not visited and it's leaf, visit
            if(!visited[tree[ancestor].children[i]] &&
               tree[tree[ancestor].children[i]].is_leaf) {
                visited[tree[ancestor].children[i]] = true;

                order_list.emplace_back(
                    tree[ancestor].children[i],
                    tree[tree[ancestor].children[i]].length + branch);
                branch = 0;
            }
        }

        // if any children is inode (not visited on foor lop above) go down that
        // branch
        if(any_of(tree[ancestor].children.begin(),
                  tree[ancestor].children.end(),
                  [&visited](int c) { return !visited[c]; })) {
            // look for !visited and inode within ancestor's children
            auto node =
                find_if(begin(tree[ancestor].children),
                        end(tree[ancestor].children), [&visited, tree](int i) {
                            return !visited[i] && !tree[i].is_leaf;
                        });

            // if so, go down that branch
            if(node != end(tree[ancestor].children)) {
                ancestor = *node;
                visited[*node] = true;
                branch += tree[ancestor].length;
            }

        } else {  // else mark ancestor as visited & move to its parent
            visited[ancestor] = true;
            branch += tree[ancestor].length;
            ancestor = tree[ancestor].parent;
        }
    }

    return EXIT_SUCCESS;
}

TEST_CASE("[tree.cc] aln_order") {
    tree_t tree;
    // cppcheck-suppress unusedVariable
    std::vector<std::pair<int, float>> order_list;
    // tree: "(B_b:6.0,(A-a:5.0,C/c:3.0,E.e:4.0)Ancestor:5.0,D%:11.0);"

    tree.emplace_back("", 0, false, 0);
    tree.emplace_back("B_b", 6, true, 0);
    tree.emplace_back("Ancestor", 5, false, 0);
    tree.emplace_back("A-a", 5, true, 2);
    tree.emplace_back("C/c", 3, true, 2);
    tree.emplace_back("E.e", 4, true, 2);
    tree.emplace_back("D%", 11, true, 0);

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
bool find_seq(const std::string& name, fasta_t& f, std::string& seq) {
    seq.clear();

    for(std::size_t i = 0; i < f.seq_names.size(); i++) {
        if(f.seq_names[i].compare(name) == 0) {
            seq = f.seq_data[i];
        }
    }

    return !seq.empty();
}

TEST_CASE("[tree.cc] find_seq") {
    // cppcheck-suppress unusedVariable
    std::string sequence;
    fasta_t fasta("", {"A", "B", "C"}, {"ACGT", "CGTA", "GTAC"});

    REQUIRE(!find_seq("Z", fasta,
                      sequence));  // fails, Z is not found -> seq is empty
    REQUIRE(find_seq("A", fasta, sequence));
    CHECK(sequence.compare("ACGT") == 0);
    REQUIRE(find_seq("B", fasta, sequence));
    CHECK(sequence.compare("CGTA") == 0);
    REQUIRE(find_seq("C", fasta, sequence));
    CHECK(sequence.compare("GTAC") == 0);
}

/* Find node in tree given its name */
bool find_node(tree_t& tree, const std::string& name, std::size_t& index) {
    auto it = find_if(begin(tree), end(tree), [name](const node_t& node) {
        return node.label == name;
    });
    index = it - begin(tree);

    return it != end(tree);
}

TEST_CASE("[tree.cc] find_node") {
    tree_t tree;
    // NOLINTNEXTLINE(clang-diagnostic-unused-variable)
    std::size_t index{0};
    // tree: "(B_b:6.0,(A-a:5.0,C/c:3.0,E.e:4.0)Ancestor:5.0,D%:11.0);"

    tree.emplace_back("", 0, false, 0);
    tree.emplace_back("B_b", 6, true, 0);
    tree.emplace_back("Ancestor", 5, false, 0);
    tree.emplace_back("A-a", 5, true, 2);
    tree.emplace_back("C/c", 3, true, 2);
    tree.emplace_back("E.e", 4, true, 2);
    tree.emplace_back("D%", 11, true, 0);

    REQUIRE(find_node(tree, tree[3].label, index));
    CHECK(index == 3);
    REQUIRE(find_node(tree, "C/c", index));
    CHECK(index == 4);
    REQUIRE(find_node(tree, "D%", index));
    CHECK(index == 6);
    REQUIRE(!find_node(tree, "Z", index));
}

/* Re-root tree given an outgroup (leaf node) */
bool reroot(tree_t& tree, const std::string& outgroup) {
    std::size_t ref{0};

    // find outgroup node
    if(!find_node(tree, outgroup, ref)) {
        throw std::invalid_argument(
            "Outgroup label could not be found, re-root failed.");
    }

    // find list of ancestors from newroot to current root
    std::vector<std::size_t> ancestors;
    std::size_t newroot = tree[ref].parent;
    std::size_t node = newroot;

    // while not current root (current root has itself as parent)
    while(tree[node].parent != node) {
        ancestors.push_back(node);
        node = tree[node].parent;
    }
    // add current root
    ancestors.push_back(node);

    // for each inode in ancestors switch the order of parent -> descendant
    for(auto i = ancestors.size() - 1; i > 0; i--) {
        tree[ancestors[i]].parent = ancestors[i - 1];
        tree[ancestors[i]].length = tree[ancestors[i - 1]].length;
    }

    // set up newroot's parent as itself with length 0
    tree[newroot].parent = newroot;
    tree[newroot].length = 0;

    return true;
}

TEST_CASE("[tree.cc] reroot") {
    tree_t tree;

    SUBCASE("One node change") {
        // tree: "(B_b:6.0,(A-a:5.0,C/c:3.0,E.e:4.0)Ancestor:5.0,D%:11.0);"
        tree.emplace_back("", 0, false, 0);
        tree.emplace_back("B_b", 6, true, 0);
        tree.emplace_back("Ancestor", 5, false, 0);
        tree.emplace_back("A-a", 5, true, 2);
        tree.emplace_back("C/c", 3, true, 2);
        tree.emplace_back("E.e", 4, true, 2);
        tree.emplace_back("D%", 11, true, 0);

        REQUIRE(reroot(tree, "A-a"));

        CHECK(tree[0].length == 5);
        CHECK(tree[0].parent == 2);
        CHECK(tree[1].length == 6);
        CHECK(tree[1].parent == 0);
        CHECK(tree[2].length == 0);
        CHECK(tree[2].parent == 2);
        CHECK(tree[3].length == 5);
        CHECK(tree[3].parent == 2);
        CHECK(tree[4].length == 3);
        CHECK(tree[4].parent == 2);
        CHECK(tree[5].length == 4);
        CHECK(tree[5].parent == 2);
        CHECK(tree[6].length == 11);
        CHECK(tree[6].parent == 0);
    }

    SUBCASE("Several node changes") {
        // tree:
        // ((raccoon:19.19959,bear:6.80041):0.84600, ((sea_lion:11.99700,
        // seal:12.00300):7.52973, ((monkey:100.85930, cat:47.14069):20.59201,
        // weasel:18.87953):2.09460):3.87382,dog:25.46154);
        tree.emplace_back("", 0, false, 0);
        tree.emplace_back("", 0.8, false, 0);
        tree.emplace_back("racoon", 19.2, true, 1);
        tree.emplace_back("bear", 6.8, true, 1);
        tree.emplace_back("", 3.9, false, 0);
        tree.emplace_back("", 7.5, false, 4);
        tree.emplace_back("sea_lion", 12, true, 5);
        tree.emplace_back("seal", 12, true, 5);
        tree.emplace_back("", 2.1, false, 4);
        tree.emplace_back("", 20.6, false, 8);
        tree.emplace_back("monkey", 100.9, true, 9);
        tree.emplace_back("cat", 47.1, true, 9);
        tree.emplace_back("weasel", 18.9, true, 8);
        tree.emplace_back("dog", 25.5, true, 0);

        REQUIRE(reroot(tree, "cat"));

        CHECK(tree[0].parent == 4);
        CHECK(tree[0].length == doctest::Approx(3.9));
        CHECK(tree[4].parent == 8);
        CHECK(tree[4].length == doctest::Approx(2.1));
        CHECK(tree[8].parent == 9);
        CHECK(tree[8].length == doctest::Approx(20.6));
        CHECK(tree[9].parent == 9);
        CHECK(tree[9].length == doctest::Approx(0));
    }
}

/* Find distance from REF to node. Tree is assumed to be rerooted. */
float distance_ref(const tree_t& tree, std::size_t ref, std::size_t node) {
    float distance = 0;

    // distance from node to root
    while(tree[node].parent != node) {
        distance += tree[node].length;
        node = tree[node].parent;
    }

    // distance from root to ref
    distance += tree[ref].length;

    return distance;
}

TEST_CASE("[tree.cc] distance_ref") {
    tree_t tree;
    // tree:
    // ((raccoon:19.19959,bear:6.80041):0.84600, ((sea_lion:11.99700,
    // seal:12.00300):7.52973, ((monkey:100.85930, cat:47.14069):20.59201,
    // weasel:18.87953):2.09460):3.87382,dog:25.46154);
    tree.emplace_back("", 0, false, 0);
    tree.emplace_back("", 0.8, false, 0);
    tree.emplace_back("racoon", 19.2, true, 1);
    tree.emplace_back("bear", 6.8, true, 1);
    tree.emplace_back("", 3.9, false, 0);
    tree.emplace_back("", 7.5, false, 4);
    tree.emplace_back("sea_lion", 12, true, 5);
    tree.emplace_back("seal", 12, true, 5);
    tree.emplace_back("", 2.1, false, 4);
    tree.emplace_back("", 20.6, false, 8);
    tree.emplace_back("monkey", 100.9, true, 9);
    tree.emplace_back("cat", 47.1, true, 9);
    tree.emplace_back("weasel", 18.9, true, 8);
    tree.emplace_back("dog", 25.5, true, 0);

    CHECK(distance_ref(tree, 13, 2) == doctest::Approx(45.5));
    CHECK(distance_ref(tree, 13, 6) == doctest::Approx(48.9));
    CHECK(distance_ref(tree, 13, 12) == doctest::Approx(50.4));
    CHECK(distance_ref(tree, 13, 11) == doctest::Approx(99.2));
}
