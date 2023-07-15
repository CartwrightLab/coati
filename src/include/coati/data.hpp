/*
# Copyright (c) 2023 Reed A. Cartwright <racartwright@gmail.com>
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

#ifndef COATI_DATA_HPP
#define COATI_DATA_HPP

#include <fst/arc.h>
#include <fst/vector-fst.h>

#include <filesystem>
#include <string>
#include <vector>

#include "matrix.hpp"

namespace coati {

using VectorFstStdArc = fst::VectorFst<fst::StdArc>;

/**
 * @brief Store sequence information.
 *
 */
class data_t {
   public:
    std::filesystem::path path;        /*!< path to input file */
    std::vector<std::string> names;    /*!< names of fasta sequences */
    std::vector<std::string> seqs;     /*!< fasta sequences */
    float_t score{0.f};                /*!< alignment score */
    std::vector<VectorFstStdArc> fsts; /*!< sequences as FSTs */
    std::vector<std::string> stops;

    data_t() = default;
    explicit data_t(std::filesystem::path p, std::vector<std::string> n = {},
                    std::vector<std::string> s = {}, float_t w = 0.f,
                    std::vector<VectorFstStdArc> f = {}, std::string c = {})
        : path{std::move(p)},
          names{std::move(n)},
          seqs{std::move(s)},
          score{w},
          fsts{std::move(f)},
          stops{std::move(c)} {}

    /** \brief Return number of names/sequences */
    [[nodiscard]] size_t size() const {
        if(names.size() != seqs.size()) {
            throw std::invalid_argument(
                "Different number of sequences and names.");
        }
        return names.size();
    }

    /** \brief Return length of sequence on position index */
    [[nodiscard]] size_t len(size_t index) const {
        return seqs[index].length();
    }
};

}  // namespace coati

#endif
