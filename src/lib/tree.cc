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

#include <boost/algorithm/string/erase.hpp>
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
// NOLINTNEXTLINE(cert-err58-cpp)
x3::rule<class label, std::string> const label = "label";
// NOLINTNEXTLINE(cert-err58-cpp)
x3::rule<class ilabel, std::string> const ilabel = "ilabel";
// NOLINTNEXTLINE(cert-err58-cpp)
x3::rule<class length, float> const length = "length";
// NOLINTNEXTLINE(cert-err58-cpp)
x3::rule<class leaf, coati::tree::tree_t> const leaf = "leaf";
// NOLINTNEXTLINE(cert-err58-cpp)
x3::rule<class inode, coati::tree::tree_t> const inode = "inode";
// NOLINTNEXTLINE(cert-err58-cpp)
x3::rule<class node, coati::tree::tree_t> const node = "node";
// NOLINTNEXTLINE(cert-err58-cpp)
x3::rule<class tree, coati::tree::tree_t> const tree = "tree";

// semantic actions
auto const make_leaf = [](auto& ctx) {
    auto label = at_c<0>(_attr(ctx));  // get label
    auto len = at_c<1>(_attr(ctx));    // get br length
    _val(ctx) = coati::tree::tree_t(
        1, {label, len, true});  // set to tree_t with 1 leaf node
};

auto const make_inode = [](auto& ctx) {
    auto label = at_c<1>(_attr(ctx));  // get label
    auto len = at_c<2>(_attr(ctx));    // get br length
    _val(ctx) =
        coati::tree::tree_t(1, {label, len});  // set to tree_t with 1 (i)node
    auto v = at_c<0>(_attr(ctx));              // get node list (tree_t)
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

namespace coati::tree {

/**
 * @brief Read tree in newick format file.
 *
 * @param[in] tree_file std::string path to newick file.
 *
 * @retval std::string content of newick file.
 */
std::string read_newick(const std::string& tree_file) {
    std::ifstream input(tree_file);  // open input stream
    if(!input.good()) {
        throw std::invalid_argument("Error opening " + tree_file + ".");
    }

    // read newick tree file
    std::string content((std::istreambuf_iterator<char>(input)),
                        std::istreambuf_iterator<char>());

    if(content.length() == 0) {  // Check file isn't empty
        throw std::invalid_argument("Reading tree failed, file is empty!");
    }

    return content;
}

/// private
// GCOVR_EXCL_START
TEST_CASE("read_newick") {
    SUBCASE("default") {
        std::ofstream outfile;
        outfile.open("tree.newick");
        REQUIRE(outfile);
        outfile << "(B_b:6.0,(A-a:5.0,C/c:3.0,E.e:4.0)Ancestor:5.0);\n";
        outfile.close();

        CHECK_EQ(read_newick("tree.newick"),
                 "(B_b:6.0,(A-a:5.0,C/c:3.0,E.e:4.0)Ancestor:5.0);\n");
        REQUIRE(std::filesystem::remove("tree.newick"));
    }
    SUBCASE("empty file - fail") {
        std::ofstream outfile;
        outfile.open("tree.newick");
        REQUIRE(outfile);
        outfile << "";
        outfile.close();

        CHECK_THROWS_AS(read_newick("tree.newick"), std::invalid_argument);
        REQUIRE(std::filesystem::remove("tree.newick"));
    }
    SUBCASE("file not found - fail") {
        CHECK_THROWS_AS(read_newick("tree.newick"), std::invalid_argument);
    }
}
// GCOVR_EXCL_STOP

/**
 * @brief Parse newick format tree.
 *
 * @details NOTE: quotation marks on labels not supported.
 *
 * @param[in] content std::string tree in newick format.
 *
 * @retval coati::tree_t parsed tree.
 */
tree_t parse_newick(std::string& content) {
    // remove tabs \t ,new lines \n, and spaces
    boost::algorithm::erase_all(content, "\t");
    boost::algorithm::erase_all(content, "\n");
    boost::algorithm::erase_all(content, " ");

    auto it = content.begin();
    auto end = content.end();

    tree_t guide_tree;
    // NOLINTNEXTLINE(misc-no-recursion)
    bool result = boost::spirit::x3::parse(it, end, newick::tree, guide_tree);

    if(!(result && it == end)) {
        throw std::runtime_error("Parsing content of newick tree failed.");
    }

    return guide_tree;
}

/// @private
// GCOVR_EXCL_START
TEST_CASE("parse_newick") {
    SUBCASE("default") {
        std::string stree{
            "(B_b:6.0,(A-a:5.0,C/c:3.0,E.e:4.0)Ancestor:5.0,D%:11.0);"};
        tree_t tree = parse_newick(stree);
        REQUIRE_EQ(tree.size(), 7);

        CHECK_EQ(tree[0].length, 0);
        CHECK_EQ(tree[0].is_leaf, false);
        CHECK_EQ(tree[0].parent, 0);

        CHECK_EQ(tree[1].label.compare("B_b"), 0);
        CHECK_EQ(tree[1].length, 6);
        CHECK_EQ(tree[1].is_leaf, true);
        CHECK_EQ(tree[1].parent, 0);

        CHECK_EQ(tree[2].label.compare("Ancestor"), 0);
        CHECK_EQ(tree[2].length, 5);
        CHECK_EQ(tree[2].is_leaf, false);
        CHECK_EQ(tree[2].parent, 0);

        CHECK_EQ(tree[3].label.compare("A-a"), 0);
        CHECK_EQ(tree[3].length, 5);
        CHECK_EQ(tree[3].is_leaf, true);
        CHECK_EQ(tree[3].parent, 2);

        CHECK_EQ(tree[4].label.compare("C/c"), 0);
        CHECK_EQ(tree[4].length, 3);
        CHECK_EQ(tree[4].is_leaf, true);
        CHECK_EQ(tree[4].parent, 2);

        CHECK_EQ(tree[5].label.compare("E.e"), 0);
        CHECK_EQ(tree[5].length, 4);
        CHECK_EQ(tree[5].is_leaf, true);
        CHECK_EQ(tree[5].parent, 2);

        CHECK_EQ(tree[6].label.compare("D%"), 0);
        CHECK_EQ(tree[6].length, 11);
        CHECK_EQ(tree[6].is_leaf, true);
        CHECK_EQ(tree[6].parent, 0);
    }
    SUBCASE("empty tree - fail") {
        std::string tree{""};
        CHECK_THROWS_AS(parse_newick(tree), std::runtime_error);
    }
}
// GCOVR_EXCL_STOP

/**
 * @brief Find sequence in data_t given its name.
 *
 * @param[in] name std::string sequence name.
 * @param[in] f coati::data_t names and sequences.
 *
 * @retval std::string sequence content.
 */
std::string find_seq(const std::string_view name, const coati::data_t& f) {
    const auto seq = std::find(f.names.cbegin(), f.names.cend(), name);

    if(seq == f.names.cend()) {
        throw std::invalid_argument("Sequence " + std::string(name) +
                                    " not found.");
    }

    return f.seqs[std::distance(f.names.cbegin(), seq)];
}

/// @private
// GCOVR_EXCL_START
TEST_CASE("find_seq") {
    // cppcheck-suppress unusedVariable
    std::string sequence;
    coati::data_t fasta("", {"A", "B", "C"}, {"ACGT", "CGTA", "GTAC"});

    sequence = find_seq("A", fasta);
    CHECK_EQ(sequence, "ACGT");
    sequence = find_seq("B", fasta);
    CHECK_EQ(sequence, "CGTA");
    sequence = find_seq("C", fasta);
    CHECK_EQ(sequence, "GTAC");
    // fails, Z is not found -> seq is empty
    CHECK_THROWS_AS(find_seq("Z", fasta), std::invalid_argument);
}
// GCOVR_EXCL_STOP

/**
 * @brief Find position of node in tree given its name.
 *
 * @param[in] tree coati::tree:tree_t phylogenetic tree.
 * @param[in] name std::string node name.
 *
 * @retval std::size_t index of node in tree.
 */
size_t find_node(const tree_t& tree, const std::string_view name) {
    auto it = find_if(begin(tree), end(tree), [name](const node_t& node) {
        return node.label == name;
    });

    if(it == end(tree)) {
        throw std::invalid_argument("Node " + std::string(name) +
                                    " not found.");
    }

    return std::distance(begin(tree), it);
}

/// @private
// GCOVR_EXCL_START
TEST_CASE("find_node") {
    tree_t tree;
    // tree: "(B_b:6.0,(A-a:5.0,C/c:3.0,E.e:4.0)Ancestor:5.0,D%:11.0);"

    tree.emplace_back("", 0, false, 0);
    tree.emplace_back("B_b", 6, true, 0);
    tree.emplace_back("Ancestor", 5, false, 0);
    tree.emplace_back("A-a", 5, true, 2);
    tree.emplace_back("C/c", 3, true, 2);
    tree.emplace_back("E.e", 4, true, 2);
    tree.emplace_back("D%", 11, true, 0);

    CHECK_EQ(find_node(tree, tree[3].label), 3);
    CHECK_EQ(find_node(tree, "C/c"), 4);
    CHECK_EQ(find_node(tree, "D%"), 6);
    REQUIRE_THROWS_AS(find_node(tree, "Z"), std::invalid_argument);
}
// GCOVR_EXCL_STOP

/**
 * @brief Re-root tree.
 *
 * @details Given the name of a leaf node, re-root the tree to set the node as
 * the outgroup.
 *
 * @param[in,out] tree coati::tree:tree_t phylogenetic tree re-rooted.
 * @param[in] nroot_name std::string node to be the new outgroup.
 */
void reroot(tree_t& tree, const std::string_view nroot_name) {
    // find new root node
    std::size_t ref = find_node(tree, nroot_name);

    // find list of ancestors from newroot to current root
    std::vector<std::size_t> ancestors;
    std::size_t newroot = tree[ref].parent;
    std::size_t node = newroot;

    // look for current root by going up the tree
    // current root has itself as parent
    while(tree[node].parent != node) {
        ancestors.push_back(node);
        node = tree[node].parent;
    }
    // add current root
    ancestors.push_back(node);

    // for each inode in ancestors swap the order of parent -> descendant
    for(auto i = ancestors.size() - 1; i > 0; i--) {
        tree[ancestors[i]].parent = ancestors[i - 1];
        tree[ancestors[i]].length = tree[ancestors[i - 1]].length;
    }

    // set up newroot's parent as itself with length 0
    tree[newroot].parent = newroot;
    tree[newroot].length = 0;
}

/// @private
// GCOVR_EXCL_START
TEST_CASE("reroot") {
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

        reroot(tree, "A-a");

        CHECK_EQ(tree[0].length, 5);
        CHECK_EQ(tree[0].parent, 2);
        CHECK_EQ(tree[1].length, 6);
        CHECK_EQ(tree[1].parent, 0);
        CHECK_EQ(tree[2].length, 0);
        CHECK_EQ(tree[2].parent, 2);
        CHECK_EQ(tree[3].length, 5);
        CHECK_EQ(tree[3].parent, 2);
        CHECK_EQ(tree[4].length, 3);
        CHECK_EQ(tree[4].parent, 2);
        CHECK_EQ(tree[5].length, 4);
        CHECK_EQ(tree[5].parent, 2);
        CHECK_EQ(tree[6].length, 11);
        CHECK_EQ(tree[6].parent, 0);
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

        reroot(tree, "cat");

        CHECK_EQ(tree[0].parent, 4);
        CHECK_EQ(tree[0].length, doctest::Approx(3.9));
        CHECK_EQ(tree[4].parent, 8);
        CHECK_EQ(tree[4].length, doctest::Approx(2.1));
        CHECK_EQ(tree[8].parent, 9);
        CHECK_EQ(tree[8].length, doctest::Approx(20.6));
        CHECK_EQ(tree[9].parent, 9);
        CHECK_EQ(tree[9].length, doctest::Approx(0));
    }
}
// GCOVR_EXCL_STOP

/**
 * @brief Find distance from reference to node.
 *
 * @details Note: tree is assumed to be re-rooted. Therefore distance is node to
 * root + root to ref (outgroup).
 *
 * @param[in] tree coati::tree::tree_t phylogenetic tree.
 * @param[in] ref std::size_t index of reference sequence in tree.
 * @param[in] node std::size_t index of node in tree.
 *
 * @retval float distance from reference sequence to node.
 */
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

/// @private
// GCOVR_EXCL_START
TEST_CASE("distance_ref") {
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
// GCOVR_EXCL_STOP
}  // namespace coati::tree
