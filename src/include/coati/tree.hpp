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

#ifndef TREE_HPP
#define TREE_HPP

#include <algorithm>
#include <boost/spirit/home/x3.hpp>
#include <string_view>

#include "utils.hpp"

namespace coati::tree {

/**
 * @brief Stores information about a tree node, used for msa.
 *
 */
struct node_t {
    std::string label; /*!< node name */
    float length; /*!< branch length connecting node to most recent ancestor */
    bool is_leaf; /*!< true if node is leaf, false otherwise */
    size_t parent{0}; /*!< node index of parent (ancestor) in tree */
    std::vector<size_t> children; /*!< children of node */

    node_t(std::string name, float len, bool leaf = false, size_t ancestor = 0)
        : label{std::move(name)},
          length{len},
          is_leaf{leaf},
          parent{ancestor} {}
};

using tree_t = std::vector<node_t>;

// Read tree in newick format file.
std::string read_newick(const std::string& tree_file);
// Parse newick format tree.
tree_t parse_newick(std::string& content);
// Find sequence in data_t given its name.
std::string find_seq(const std::string_view name, const coati::data_t& f);
// Find position of node in tree given its name.
size_t find_node(const tree_t& tree, const std::string_view name);
// Re-root tree.
void reroot(tree_t& tree, const std::string_view label);
// Find distance from reference to node.
float distance_ref(const tree_t& tree, size_t ref, size_t node);
}  // namespace coati::tree
#endif
