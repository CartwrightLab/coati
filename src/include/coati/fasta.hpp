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

#ifndef FASTA_HPP
#define FASTA_HPP

#include <fst/fstlib.h>

#include <filesystem>
#include <string>
#include <vector>

using VectorFstStdArc = fst::VectorFst<fst::StdArc>;

namespace coati {
struct fasta_t {
    std::filesystem::path path;
    std::vector<std::string> names, seqs;

    fasta_t() = default;
    explicit fasta_t(std::filesystem::path p, std::vector<std::string> n = {},
                     std::vector<std::string> s = {})
        : path{std::move(p)}, names{std::move(n)}, seqs{std::move(s)} {}
    size_t size() {
        if(names.size() != seqs.size()) {
            throw std::invalid_argument(
                "Different number of sequences and names in fasta file.");
        }
        return names.size();
    }
    [[nodiscard]] size_t size() const {
        if(names.size() != seqs.size()) {
            throw std::invalid_argument(
                "Different number of sequences and names in fasta file.");
        }
        return names.size();
    }
};

fasta_t read_fasta(const std::string& f_path,
                   std::vector<VectorFstStdArc>& fsts);
fasta_t read_fasta(const std::string& f_path);
bool write_fasta(const VectorFstStdArc& aln, fasta_t& fasta);
bool write_fasta(const fasta_t& fasta);

}  // namespace coati

#include "mutation_fst.hpp"

#endif
