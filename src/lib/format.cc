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
int format_sequences(coati::format_t& format, coati::alignment_t& aln) {
    if(format.preserve_phase) {  // padd gaps to preserve phase
        if(format.padding == "-") {
            throw std::invalid_argument("Invalid padding character " +
                                        format.padding + " .");
        }
        auto pos = aln.seq(0).find('-');
        while(pos != std::string::npos) {
            size_t len{0};
            while(aln.seq(0)[pos] == '-') {
                pos++;
                len++;
            }
            len = len % 3;
            for(size_t i = 0; i < aln.data.size(); i++) {
                switch(len) {
                case 1:
                    aln.seq(i).insert(pos, format.padding, 0, len);
                case 2:
                    aln.seq(i).insert(pos, format.padding, 0, len);
                }
            }
            pos = aln.seq(0).find('-', pos++);
        }
    }

    // if sequences to extract are specified do so
    if(format.seqs.size() > 0 || format.pos.size() > 0) {
        for(auto index : format.pos) {
            format.seqs.emplace_back(aln.data.names[index - 1]);
        }

        for(int i = aln.data.names.size() - 1; i >= 0; --i) {
            // if sequence name is not found on extraction list, remove
            if(std::find(begin(format.seqs), end(format.seqs),
                         aln.data.names[i]) == end(format.seqs)) {
                aln.data.names.erase(aln.data.names.begin() + i);
                aln.data.seqs.erase(aln.data.seqs.begin() + i);
            }
        }

        // if all sequences are removed, throw an error
        if(aln.data.names.empty()) {
            throw std::invalid_argument("Sequences not found.");
        }
    }

    // output formatted sequences
    return utils::write_output(aln.data) ? 0 : 1;
}

/// @private
// GCOVR_EXCL_START
TEST_CASE("format_sequences") {
    coati::args_t args;
    coati::data_t expected;

    // NOLINTNEXTLINE(misc-unused-parameters)
    auto test_format_fasta = [](args_t args, data_t expected) {
        args.aln.data.out_file = {"test-format-seqs.fa", ".fa"};
        if(std::filesystem::exists(args.aln.data.out_file.path)) {
            std::filesystem::remove(args.aln.data.out_file.path);
        }

        REQUIRE(coati::format_sequences(args.format, args.aln) == 0);

        std::ifstream infile(args.aln.data.out_file.path);
        std::string s1, s2;

        for(size_t i = 0; i < expected.size(); ++i) {
            infile >> s1 >> s2;
            CHECK(s1 == expected.names[i]);
            CHECK(s2 == expected.seqs[i]);
        }

        CHECK(std::filesystem::remove(args.aln.data.out_file.path));
    };

    // NOLINTNEXTLINE(misc-unused-parameters)
    auto test_format_phylip = [](args_t args, data_t expected) {
        args.aln.data.out_file = {"test-format-seqs.phy", ".phy"};
        if(std::filesystem::exists(args.aln.data.out_file.path)) {
            std::filesystem::remove(args.aln.data.out_file.path);
        }

        REQUIRE(coati::format_sequences(args.format, args.aln) == 0);

        std::ifstream infile(args.aln.data.out_file.path);
        std::string s1, s2;

        infile >> s1 >> s2;
        CHECK(std::stoul(s1) == expected.size());
        CHECK(std::stoul(s2) == args.aln.seq(0).size());

        for(size_t i = 0; i < expected.size(); ++i) {
            infile >> s1 >> s2;
            CHECK(s1 == expected.names[i]);
            CHECK(s2 == expected.seqs[i]);
        }

        CHECK(std::filesystem::remove(args.aln.data.out_file.path));
    };

    SUBCASE("fasta-multiple-3") {
        args.format.preserve_phase = true;
        args.aln.data =
            coati::data_t("", {"Ancestor", "Descendant"}, {"AG---T", "ACCCGT"});
        expected = coati::data_t("", {">Ancestor", ">Descendant"},
                                 {"AG---T", "ACCCGT"});
        test_format_fasta(args, expected);
    }
    SUBCASE("phylip-3seqs-no-preserve-case") {
        args.aln.data = coati::data_t(
            "", {"Ancestor", "Descend-1", "Descend-2"}, {"AGT", "AGT", "AG-"});
        expected = coati::data_t("", {"Ancestor", "Descend-1", "Descend-2"},
                                 {"AGT", "AGT", "AG-"});
        test_format_phylip(args, expected);
    }
    SUBCASE("fasta-gap-len2") {
        args.format.preserve_phase = true;
        args.aln.data =
            coati::data_t("", {"Ancestor", "Descendant"}, {"A--GT", "ACCGT"});
        expected = coati::data_t("", {">Ancestor", ">Descendant"},
                                 {"A--?GT", "ACC?GT"});
        test_format_fasta(args, expected);
    }
    SUBCASE("fasta-gap-len11") {
        args.format.preserve_phase = true;
        args.aln.data = coati::data_t("", {"Ancestor", "Descendant"},
                                      {"A-----------GT", "ACCCCCCCCCCCGT"});
        expected = coati::data_t("", {">Ancestor", ">Descendant"},
                                 {"A-----------?GT", "ACCCCCCCCCCC?GT"});
        test_format_fasta(args, expected);
    }
    SUBCASE("phylip-3seqs-gap-len1") {
        args.format.preserve_phase = true;
        args.format.padding = "X";
        args.aln.data =
            coati::data_t("", {"Ancestor", "Descend-1", "Descend-2"},
                          {"A-GT", "ACGT", "ACG-"});
        expected = coati::data_t("", {"Ancestor", "Descend-1", "Descend-2"},
                                 {"A-XXGT", "ACXXGT", "ACXXG-"});
        test_format_phylip(args, expected);
    }
    SUBCASE("fasta-extract-sequences-name") {
        args.format.preserve_phase = true;
        args.aln.data = coati::data_t("", {"Ancestor", "Descendant"},
                                      {"A-----------GT", "ACCCCCCCCCCCGT"});
        args.format.seqs = {"Ancestor"};
        expected = coati::data_t("", {">Ancestor"}, {"A-----------?GT"});
        test_format_fasta(args, expected);
    }
    SUBCASE("fasta-extract-sequences-position") {
        args.format.preserve_phase = true;
        args.aln.data =
            coati::data_t("", {"Ancestor", "Descendant1", "Descendant2"},
                          {"A--GT", "ACCGT", "A-CGT"});
        args.format.pos = {1, 3};
        expected = coati::data_t("", {">Ancestor", ">Descendant2"},
                                 {"A--?GT", "A-C?GT"});
        test_format_fasta(args, expected);
    }
    SUBCASE("phylip-extract-sequences-name-position") {
        args.format.preserve_phase = true;
        args.aln.data =
            coati::data_t("", {"Ancestor", "Descend-1", "Descend-2"},
                          {"A--GT", "ACCGT", "A-CGT"});
        args.format.seqs = {"Ancestor"};
        args.format.pos = {1, 3};
        expected =
            coati::data_t("", {"Ancestor", "Descend-2"}, {"A--?GT", "A-C?GT"});
        test_format_phylip(args, expected);
    }
    SUBCASE("phylip-extract-sequences-name-position2") {
        args.format.preserve_phase = true;
        args.format.padding = "$";
        args.aln.data =
            coati::data_t("", {"Ancestor", "Descend-1", "Descend-2"},
                          {"A-GT", "ACGT", "AC-T"});
        args.format.seqs = {"Ancestor"};
        args.format.pos = {2};
        expected =
            coati::data_t("", {"Ancestor", "Descend-1"}, {"A-$$GT", "AC$$GT"});
        test_format_phylip(args, expected);
    }
    SUBCASE("invalid-padding") {
        args.format.preserve_phase = true;
        args.format.padding = "-";
        REQUIRE_THROWS_AS(coati::format_sequences(args.format, args.aln),
                          std::invalid_argument);
    }
    SUBCASE("sequences-not-found") {
        args.aln.data =
            coati::data_t("", {"Ancestor", "Descendant"}, {"AAA", "AAA"});
        args.format.seqs = {"coati"};
        REQUIRE_THROWS_AS(coati::format_sequences(args.format, args.aln),
                          std::invalid_argument);
    }
}
// GCOVR_EXCL_STOP
}  // namespace coati
