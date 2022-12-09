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
 * @brief Format and manage sequences.
 *
 * @details Format input sequences to user specified format (FASTA, PHYLIP,
 * JSON), add padding to preserve ancestor phase, extract and/or reorder
 * sequences.
 *
 * @param[in,out] format coati::format_t format arguments.
 * @param[in,out] aln coati::alignment_t alignment arguments.
 *
 * @retval int EXIT_SUCCESS (0) if no errors.
 */
int format_sequences(coati::format_t& format, coati::alignment_t& aln) {
    // if sequences to extract are specified do so
    if(format.names.size() > 0 || format.pos.size() > 0) {
        extract_seqs(format, aln.data);
    }

    if(format.preserve_phase) {      // padd gaps to preserve phase
        if(format.padding == "-") {  // padding char cannot be same as gap char
            throw std::invalid_argument("Invalid padding character " +
                                        format.padding + " .");
        }
        auto pos = aln.seq(0).find('-');  // find first gap (insertion)
        while(pos != std::string::npos) {
            size_t len{0};
            while(aln.seq(0)[pos] == '-') {  // go to end of gap
                pos++;
                len++;
            }
            len = len % 3;
            // add 1 or 2 padding chars so that next codon starts in frame
            for(size_t i = 0; i < aln.data.size(); i++) {
                switch(len) {
                case 1:
                    aln.seq(i).insert(pos, format.padding, 0, len);
                case 2:
                    aln.seq(i).insert(pos, format.padding, 0, len);
                }
            }
            pos = aln.seq(0).find('-', pos++);  // find next gap (insertion)
        }
    }

    // output formatted sequences
    coati::io::write_output(aln);
    return EXIT_SUCCESS;
}

/**
 * @brief Keep only sequences specified by position or by name.
 *
 * @details Remove all sequences that are not specified and reorder them if
 * necessary. Sequence to extract can be specified by either position or name,
 * but not both simultaneously.
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

    // if specified by name
    if(format.names.size() > 0) {
        for(std::size_t i = 0; i < format.names.size(); ++i) {
            auto it =
                find(data.names.cbegin(), data.names.cend(), format.names[i]);
            if(it != data.names.cend()) {
                // add +1 so that pos are consistent with user input (1 indexed)
                format.pos.push_back(std::distance(data.names.cbegin(), it) +
                                     1);
            } else {  // sequence not found
                throw std::invalid_argument("Sequence " + format.names[i] +
                                            " not found.");
            }
        }
    }

    // if specified by position and extract reorder when specified by name
    if(format.pos.size() > 0) {
        // check that positions of sequence are valid (not out of bounds)
        const auto [min, max] =
            std::minmax_element(format.pos.cbegin(), format.pos.cend());
        if(*min == 0 || *max > data.size()) {
            throw std::invalid_argument(
                "Positions of seqs to extract are of out range");
        }
        // extract sequences and their names
        for(auto pos : format.pos) {
            names.push_back(data.names[pos - 1]);
            seqs.push_back(data.seqs[pos - 1]);
        }
    }
    data.names = names;
    data.seqs = seqs;
}
/// @private
// GCOVR_EXCL_START
TEST_CASE("extract_seqs") {
    auto test = [](coati::format_t& format, coati::data_t& data,
                   // NOLINTNEXTLINE(misc-unused-parameters)
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
    SUBCASE("Select one seq by pos") {
        coati::data_t data, expected;
        data.names = {"A", "B"};
        data.seqs = {"AAA", "CCC"};
        coati::format_t format;
        format.pos = {2};

        expected.names = {"B"};
        expected.seqs = {"CCC"};

        test(format, data, expected);
    }
    SUBCASE("Select two seqs: one found, one not - fail") {
        coati::data_t data, expected;
        data.names = {"A", "B", "C"};
        data.seqs = {"AAA", "GGG", "CCC"};
        coati::format_t format;
        format.names = {"C", "D"};

        CHECK_THROWS_AS(extract_seqs(format, data), std::invalid_argument);
    }
    SUBCASE("Names not found - fail") {
        coati::data_t data, expected;
        data.names = {"A", "B"};
        data.seqs = {"AAA", "CCC"};

        coati::format_t format;
        format.names = {"C"};

        CHECK_THROWS_AS(extract_seqs(format, data), std::invalid_argument);
    }
    SUBCASE("Positions out of range - fail") {
        coati::data_t data, expected;
        data.names = {"A", "B"};
        data.seqs = {"AAA", "CCC"};

        coati::format_t format;
        format.pos = {5};

        CHECK_THROWS_AS(extract_seqs(format, data), std::invalid_argument);
    }
    SUBCASE("Positions out of range - fail") {
        coati::data_t data, expected;
        data.names = {"A", "B"};
        data.seqs = {"AAA", "CCC"};

        coati::format_t format;
        format.pos = {0};

        CHECK_THROWS_AS(extract_seqs(format, data), std::invalid_argument);
    }
}
// GCOVR_EXCL_STOP

/// @private
// GCOVR_EXCL_START
TEST_CASE("format_sequences") {
    coati::args_t args;
    coati::data_t expected;

    // NOLINTNEXTLINE(misc-unused-parameters)
    auto test_format_fasta = [](args_t args, data_t expected) {
        args.aln.output = "test-format-seqs.fa";
        if(std::filesystem::exists(args.aln.output)) {
            std::filesystem::remove(args.aln.output);
        }

        REQUIRE_EQ(coati::format_sequences(args.format, args.aln), 0);

        std::ifstream infile(args.aln.output);
        std::string s1, s2;

        for(size_t i = 0; i < expected.size(); ++i) {
            infile >> s1 >> s2;
            CHECK(s1 == expected.names[i]);
            CHECK(s2 == expected.seqs[i]);
        }

        CHECK(std::filesystem::remove(args.aln.output));
    };

    // NOLINTNEXTLINE(misc-unused-parameters)
    auto test_format_phylip = [](args_t args, data_t expected) {
        args.aln.output = "test-format-seqs.phy";
        if(std::filesystem::exists(args.aln.output)) {
            std::filesystem::remove(args.aln.output);
        }

        REQUIRE_EQ(coati::format_sequences(args.format, args.aln), 0);

        std::ifstream infile(args.aln.output);
        std::string s1, s2;

        infile >> s1 >> s2;
        CHECK(std::stoul(s1) == expected.size());
        CHECK(std::stoul(s2) == args.aln.seq(0).size());

        for(size_t i = 0; i < expected.size(); ++i) {
            infile >> s1 >> s2;
            CHECK(s1 == expected.names[i]);
            CHECK(s2 == expected.seqs[i]);
        }

        CHECK(std::filesystem::remove(args.aln.output));
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
