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

#include <coati/format.hpp>

namespace coati {

/**
 * \brief Format sequences according to user input.
 *
 * Format input sequences to format specified by user (FASTA, PHYLIP, JSON)
 *  and add padding to preserve ancestor phase if required.
 *
 * @param[in] args coati::args_t input arguments;
 *
 */
int format_sequences(coati::args_t& args) {
    if(args.preserve_phase) {  // padd gaps to preserve phase
        if(args.padding == "-") {
            throw std::invalid_argument("Invalid padding character " +
                                        args.padding + " .");
        }
        auto pos = args.aln.data.seqs[0].find('-');
        while(pos != std::string::npos) {
            size_t len{0};
            while(args.aln.data.seqs[0][pos] == '-') {
                pos++;
                len++;
            }
            for(size_t i = 0; i < args.aln.data.size(); i++) {
                switch(len) {
                case 1:
                    args.aln.data.seqs[i].insert(pos, args.padding, 0, len);
                case 2:
                    args.aln.data.seqs[i].insert(pos, args.padding, 0, len);
                }
            }
            pos = args.aln.data.seqs[0].find('-', pos++);
        }
    }

    // output formatted sequences
    return utils::write_output(args.aln.data) ? 0 : 1;
}

/// @private
TEST_CASE("format_sequences") {
    coati::args_t args;

    SUBCASE("fasta") {
        args.preserve_phase = true;
        args.aln.data =
            coati::data_t("", {"Ancestor", "Descendant"}, {"AG---T", "ACCCGT"});
        args.aln.data.out_file = {"test-format-fasta-fasta.fa", ".fa"};
        REQUIRE(coati::format_sequences(args) == 0);
        // read file
        std::ifstream infile(args.aln.data.out_file.path);
        std::string s1, s2;

        infile >> s1 >> s2;
        CHECK(s1 == ">Ancestor");
        CHECK(s2 == "AG---T");

        infile >> s1 >> s2;
        CHECK(s1 == ">Descendant");
        CHECK(s2 == "ACCCGT");
        // delete file
        CHECK(std::filesystem::remove(args.aln.data.out_file.path));
    }
    SUBCASE("phylip") {
        args.aln.data =
            coati::data_t("", {"Ancestor", "Descendant1", "Descendant2"},
                          {"AGT", "AGT", "AG-"});
        args.aln.data.out_file = {"test-format-fasta-phylip.phy", ".phy"};
        REQUIRE(coati::format_sequences(args) == 0);
        // read file
        std::ifstream infile(args.aln.data.out_file.path);
        std::string s1, s2;

        infile >> s1 >> s2;
        CHECK(s1 == "3");
        CHECK(s2 == "6");

        infile >> s1 >> s2;
        CHECK(s1 == "Ancestor");
        CHECK(s2 == "AGT");

        infile >> s1 >> s2;
        CHECK(s1 == "Descendant1");
        CHECK(s2 == "AGT");

        infile >> s1 >> s2;
        CHECK(s1 == "Descendant2");
        CHECK(s2 == "AG-");
        // delete file
        CHECK(std::filesystem::remove(args.aln.data.out_file.path));
    }
    SUBCASE("fasta") {
        args.preserve_phase = true;
        args.aln.data =
            coati::data_t("", {"Ancestor", "Descendant"}, {"A--GT", "ACCGT"});
        args.aln.data.out_file = {{"test-format-fasta.fa", ".fa"}};
        REQUIRE(coati::format_sequences(args) == 0);
        // read file
        std::ifstream infile(args.aln.data.out_file.path);
        std::string s1, s2;

        infile >> s1 >> s2;
        CHECK(s1 == ">Ancestor");
        CHECK(s2 == "A--?GT");

        infile >> s1 >> s2;
        CHECK(s1 == ">Descendant");
        CHECK(s2 == "ACC?GT");
        // delete file
        CHECK(std::filesystem::remove(args.aln.data.out_file.path));
    }
    SUBCASE("phylip") {
        args.preserve_phase = true;
        args.padding = "X";
        args.aln.data =
            coati::data_t("", {"Ancestor", "Descendant1", "Descendant2"},
                          {"A-GT", "ACGT", "ACG-"});
        args.aln.data.out_file = {{"test-format-phylip.phy", ".phy"}};
        REQUIRE(coati::format_sequences(args) == 0);
        // read file
        std::ifstream infile(args.aln.data.out_file.path);
        std::string s1, s2;

        infile >> s1 >> s2;
        CHECK(s1 == "3");
        CHECK(s2 == "6");

        infile >> s1 >> s2;
        CHECK(s1 == "Ancestor");
        CHECK(s2 == "A-??GT");

        infile >> s1 >> s2;
        CHECK(s1 == "Descendant1");
        CHECK(s2 == "AC??GT");

        infile >> s1 >> s2;
        CHECK(s1 == "Descendant2");
        CHECK(s2 == "AC??G-");
        // delete file
        CHECK(std::filesystem::remove(args.aln.data.out_file.path));
    }
    SUBCASE("invalid padding") {
        args.preserve_phase = true;
        args.padding = "-";
        REQUIRE_THROWS_AS(coati::format_sequences(args), std::invalid_argument);
    }
}
}  // namespace coati
