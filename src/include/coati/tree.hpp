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

#ifndef TREE_HPP
#define TREE_HPP

#include <algorithm>
#include <boost/spirit/home/x3.hpp>

#include "utils.hpp"

struct node_t {
    std::string label;
    float length;
    bool is_leaf;
    size_t parent{0};
    std::vector<int> children;

    node_t(std::string name, float len, bool leaf = false, size_t ancestor = 0)
        : label{std::move(name)},
          length{len},
          is_leaf{leaf},
          parent{ancestor} {}
};

using tree_t = std::vector<node_t>;

bool read_newick(const std::string& tree_file, std::string& content);
int parse_newick(std::string content, tree_t& guide_tree);
int aln_order(tree_t& tree, std::vector<std::pair<int, float>>& order_list);
bool find_seq(const std::string& name, fasta_t& f, std::string& seq);
bool find_node(tree_t& tree, const std::string& name, int& ID);
bool reroot(tree_t& tree, const std::string& label);
float distance_ref(const tree_t& tree, size_t ref, size_t node);

#endif
