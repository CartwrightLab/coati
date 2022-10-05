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
 * @param[in,out] format coati::format_t format arguments.
 * @param[in,out] aln coati::alignment_t alignment arguments.
 *
 */
int format_sequences(coati::format_t& format, coati::alignment_t& aln) {
    // if sequences to extract are specified do so
    if(format.names.size() > 0 || format.pos.size() > 0) {
        extract_seqs(format, aln.data);
    }

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

    // output formatted sequences
    coati::io::write_output(aln.data);
    return EXIT_SUCCESS;
}

/**
 * \brief Keep only sequences specified by position or by name.
 *
 * @param[in,out] format coati::format_t arguments for running coati-format.
 * @param[in,out] data coati::data_t information about sequence data and
 * paramaters.
 */
void extract_seqs(coati::format_t& format, coati::data_t& data) {
    std::vector<std::string> names;
    std::vector<std::string> seqs;
    names.reserve(data.names.size());
    seqs.reserve(data.seqs.size());

    if(format.names.size() > 0) {
        for(std::size_t i = 0; i < format.names.size(); ++i) {
            auto it =
                find(data.names.cbegin(), data.names.cend(), format.names[i]);
            if(it != data.names.cend()) {
                format.pos.push_back(it - data.names.cbegin() + 1);
                // add +1 so that numbers are 1 indexed - consistent with user
                //  input
            }
        }
    }

    if(format.pos.size() > 0) {
        for(auto pos : format.pos) {
            names.push_back(data.names[pos - 1]);
            seqs.push_back(data.seqs[pos - 1]);
        }
    }
    // if all sequences are removed, throw an error
    if(data.names.empty()) {
        throw std::invalid_argument("Sequences to extract not found.");
    }
    data.names = names;
    data.seqs = seqs;
}
/// @private
// GCOVR_EXCL_START
TEST_CASE("extract_seqs") {
    // NOLINTNEXTLINE(mis-unused-parameters)
    auto test = [](coati::format_t& format, coati::data_t& data,
                   const coati::data_t& expected) {
        coati::extract_seqs(format, data);
        CHECK_EQ(data.names, expected.names);
        CHECK_EQ(data.seqs, expected.seqs);
    };

    SUBCASE("Flip two sequences by name") {
        coati::data_t data, expected;
        data.names = {"A", "B"};
        data.seqs = {"AAA", "CCC"};
        coati::format_t format;
        format.names = {"B", "A"};

        expected.names = {"B", "A"};
        expected.seqs = {"CCC", "AAA"};

        test(format, data, expected);
    }

    SUBCASE("Select one seq by pos") {}
    SUBCASE("Remove all explicit - fail") {}
    SUBCASE("Remove all names not found - fail") {}
}
// GCOVR_EXCL_STOP

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

        REQUIRE_EQ(coati::format_sequences(args.format, args.aln), 0);

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

        REQUIRE_EQ(coati::format_sequences(args.format, args.aln), 0);

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
        args.format.names = {"Ancestor"};
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
    SUBCASE("phylip-extract-sequences-position") {
        args.format.preserve_phase = true;
        args.aln.data =
            coati::data_t("", {"Ancestor", "Descend-1", "Descend-2"},
                          {"A--GT", "ACCGT", "A-CGT"});
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
                          {"ACGT", "A-GT", "AC-T"});
        args.format.pos = {2, 1};
        expected =
            coati::data_t("", {"Descend-1", "Ancestor"}, {"A-$$GT", "AC$$GT"});
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
        args.format.names = {"coati"};
        REQUIRE_THROWS_AS(coati::format_sequences(args.format, args.aln),
                          std::invalid_argument);
    }
}
// GCOVR_EXCL_STOP
}  // namespace coati
