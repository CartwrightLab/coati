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

#include <coati/align_marginal.hpp>
#include <filesystem>

namespace coati {

/**
 * \brief Pairwise alignment using dynamic programming and a marginal model.
 *
 * Alignment of two sequences via dynamic programming using an affine
 *  (geometrical) gap model and a marginal version of Muse and Gaut (1984)
 *  codon substitution model.
 *
 * @param[in] aln coati::alignment_t alignment information.
 */
bool marg_alignment(coati::alignment_t& aln) {
    coati::Matrixf P(64, 64), p_marg;
    std::ofstream out_w;

    // set substitution matrix according to model
    coati::utils::set_subst(aln);

    // score alignment
    if(aln.score) {
        std::cout << alignment_score(aln, aln.subst_matrix) << std::endl;
        return true;
    }

    // check that length of ref sequence is multiple of 3 and gap unit size
    size_t len_a = aln.seq(0).length();
    if((len_a % 3 != 0) || (len_a % aln.gap.len != 0)) {
        throw std::invalid_argument(
            "Length of reference sequence must be multiple of 3.");
    }
    // check that length of descendant sequence is multiple of gap unit size
    if(aln.seq(1).length() % aln.gap.len != 0) {
        throw std::invalid_argument(
            "Length of descendant sequence must be multiple of " +
            std::to_string(aln.gap.len) + ".");
    }

    // encode sequences
    auto anc = aln.seq(0);
    auto des = aln.seq(1);
    coati::utils::sequence_pair_t seq_pair =
        coati::utils::marginal_seq_encoding(anc, des);

    // dynamic programming pairwise alignment and traceback
    coati::align_pair_work_mem_t work;
    try {
        coati::viterbi_mem(work, seq_pair[0], seq_pair[1], aln);
    } catch(const std::bad_alloc& e) {
        std::cout << "Erorr: sequences to align exceed available memory."
                  << std::endl;
        return EXIT_FAILURE;
    }
    coati::traceback(work, anc, des, aln, aln.gap.len);

    if(!aln.weight_file.empty()) {  // save weight and filename
        out_w.open(aln.weight_file, std::ios::app | std::ios::out);
        out_w << aln.data.path << "," << aln.model << "," << aln.data.weight
              << std::endl;
        out_w.close();
    }

    // write alignment
    coati::utils::write_output(aln.data);
    return true;
}

/// @private
// GCOVR_EXCL_START
TEST_CASE("marg_alignment") {
    // NOLINTNEXTLINE(misc-unused-parameters)
    auto test_fasta = [](alignment_t& aln, const data_t& expected) {
        if(std::filesystem::exists(aln.data.out_file.path)) {
            std::filesystem::remove(aln.data.out_file.path);
        }
        if(std::filesystem::exists(aln.weight_file)) {
            std::filesystem::remove(aln.weight_file);
        }
        REQUIRE(marg_alignment(aln));

        std::ifstream infile(aln.data.out_file.path);
        std::string s1, s2;

        infile >> s1 >> s2;
        CHECK(s1 == expected.names[0]);
        CHECK(s2 == expected.seqs[0]);

        infile >> s1 >> s2;
        CHECK(s1 == expected.names[1]);
        CHECK(s2 == expected.seqs[1]);
        CHECK(std::filesystem::remove(aln.data.out_file.path));

        if(!aln.weight_file.empty()) {
            std::ifstream inweight(aln.weight_file);
            std::string s;  // NOLINT(clang-diagnostic-unused-variable)
            inweight >> s;
            // NOLINTNEXTLINE(clang-diagnostic-unused-variable)
            std::size_t start = s.find_last_of(',');
            CHECK(std::filesystem::remove(aln.weight_file));
            CHECK(std::stof(s.substr(start + 1)) == expected.weight);
        }
    };

    // NOLINTNEXTLINE(misc-unused-parameters)
    auto test_phylip = [](alignment_t& aln, const data_t& expected) {
        if(std::filesystem::exists(aln.data.out_file.path)) {
            std::filesystem::remove(aln.data.out_file.path);
        }
        if(std::filesystem::exists(aln.weight_file)) {
            std::filesystem::remove(aln.weight_file);
        }
        REQUIRE(marg_alignment(aln));

        std::ifstream infile(aln.data.out_file.path);
        std::string s1, s2;

        infile >> s1 >> s2;
        CHECK(s1 == "2");
        CHECK(std::stoul(s2) == aln.seq(0).size());

        infile >> s1 >> s2;
        CHECK(s1 == expected.names[0]);
        CHECK(s2 == expected.seqs[0]);

        infile >> s1 >> s2;
        CHECK(s1 == expected.names[1]);
        CHECK(s2 == expected.seqs[1]);
        CHECK(std::filesystem::remove(aln.data.out_file.path));
    };

    alignment_t aln;
    data_t expected;
    SUBCASE("Alignment - output fasta") {
        aln.data = data_t("", {"1", "2"}, {"CTCTGGATAGTG", "CTATAGTG"});
        aln.model = "m-coati";
        aln.weight_file = "score.log";
        aln.data.out_file = {{"test-marg_alignment-fasta.fasta"}, {".fasta"}};

        expected = data_t("", {">1", ">2"}, {"CTCTGGATAGTG", "CT----ATAGTG"});
        expected.weight = 1.51294f;

        test_fasta(aln, expected);
    }
    SUBCASE("Alignment - output phylip") {
        aln.data = coati::data_t("", {"1", "2"}, {"GCGACTGTT", "GCGATTGCTGTT"});
        aln.model = "m-coati";
        aln.data.out_file = {{"test-marg_alignment-phylip.phy"}, {".phy"}};

        expected = data_t("", {"1", "2"}, {"GCGA---CTGTT", "GCGATTGCTGTT"});

        test_phylip(aln, expected);
    }
    SUBCASE("Alignment 2 dels - output phylip") {
        aln.data = coati::data_t("", {"1", "2"}, {"ACGTTAAGGGGT", "ACGAAT"});
        aln.model = "m-coati";
        aln.data.out_file = {{"test-marg_alignment-phylip2.phy"}, {".phy"}};

        expected = data_t("", {"1", "2"}, {"ACGTTAAGGGGT", "ACG--AA----T"});

        test_phylip(aln, expected);
    }
    SUBCASE("Alignment with gap length multiple of 3") {
        aln.data = coati::data_t("", {"1", "2"}, {"ACGTTAAGGGGT", "ACGAAT"});
        aln.model = "m-coati";
        aln.data.out_file = {{"test-marg_alignment-no-frameshifts.fa"},
                             {".fa"}};
        aln.gap.len = 3;

        expected =
            data_t("", {">1", ">2"}, {"ACG---TTAAGGGGT", "ACGAAT---------"});

        test_fasta(aln, expected);
    }
    SUBCASE("Alignment with ambiguous nucleotides") {
        aln.data = coati::data_t("", {"1", "2"}, {"CTCTGGATAGTG", "CTATAGTR"});
        aln.model = "m-coati";
        aln.weight_file = "score.log";
        aln.data.out_file = {{"test-marg_alignment_ambiguous.fa"}, {".fa"}};

        expected = data_t("", {">1", ">2"}, {"CTCTGGATAGTG", "CT----ATAGTR"});
        expected.weight = -1.03892f;  // AVG
        // expected.weight = 1.51294f;  // BEST
        test_fasta(aln, expected);
    }
    SUBCASE("Alignment with gap length multiple of 3 - fail") {
        aln.data = coati::data_t("", {"1", "2"}, {"GCGATTGCTGT", "GCGACTGTT"});
        aln.model = "m-coati";
        aln.data.out_file = {{"test-marg_alignment-no-frameshifts-f.fasta"},
                             {".fasta"}};
        aln.gap.len = 3;

        REQUIRE_THROWS_AS(marg_alignment(aln), std::invalid_argument);
    }
    SUBCASE("Length of descendant not multiple of gap.len") {
        aln.data = coati::data_t("", {"A", "B"}, {"CTCGGA", "CTCGG"});
        aln.model = "m-coati";
        aln.gap.len = 3;
        REQUIRE_THROWS_AS(marg_alignment(aln), std::invalid_argument);
    }
    SUBCASE("Score alignment") {
        aln.data =
            coati::data_t("", {"1", "2"}, {"CTCTGGATAGTG", "CT----ATAGTG"});
        aln.model = "m-coati";
        aln.weight_file = "";
        aln.data.out_file = {{"test-marg_alignment-score.fasta"}, {".fasta"}};
        aln.score = true;

        REQUIRE(marg_alignment(aln));
    }
    SUBCASE("Score alignment - fail") {
        aln.data = coati::data_t("", {"1", "2"}, {"CTCTGGATAGTG", "CTATAGTG"});
        aln.model = "m-coati";
        aln.weight_file = "";
        aln.data.out_file = {{"test-marg_alignment-score-f.fasta"}, {".fasta"}};
        aln.score = true;

        coati::data_t result(aln.data.out_file.path);

        REQUIRE_THROWS_AS(marg_alignment(aln), std::invalid_argument);
    }
}
// GCOVR_EXCL_STOP

/**
 * \brief Score alignment using marginal model.
 *
 * @param[in] args coati::args_t input parameters.
 * @param[in] p_marg coati::Matrixf substitution matrix.
 *
 * \return alignment score (float).
 */
float alignment_score(const coati::alignment_t& aln,
                      const coati::Matrixf& p_marg) {
    std::vector<std::string> seqs = aln.data.seqs;

    // check that both sequences have equal length
    if(seqs[0].length() != seqs[1].length()) {
        throw std::invalid_argument(
            "For alignment scoring both sequences must have equal length.");
    }

    // encode desc and gap-less ref sequences for subsitution matrix access
    std::string anc{seqs[0]};
    boost::erase_all(anc, "-");
    coati::utils::sequence_pair_t seq_pair =
        coati::utils::marginal_seq_encoding(anc, seqs[1]);

    coati::semiring::tropical tropical_rig;
    // calculate log(1-g) log(1-e) log(g) log(e) log(pi)
    float_t no_gap = tropical_rig.from_linear_1mf(aln.gap.open);
    float_t gap_stop = tropical_rig.from_linear_1mf(aln.gap.extend);
    float_t gap_open = tropical_rig.from_linearf(aln.gap.open);
    float_t gap_extend = tropical_rig.from_linearf(aln.gap.extend);
    std::vector<coati::float_t> pi{aln.pi};
    std::transform(pi.cbegin(), pi.cend(), pi.begin(),
                   [](auto value) { return ::logf(value); });

    float weight{0.f};
    int state{0}, ngap{0};
    for(size_t i = 0; i < seqs[0].length(); i++) {
        switch(state) {
        case 0:  // subsitution
            if(seqs[0][i] == '-') {
                // insertion;
                weight += gap_open;
                state = 2;
                ngap++;
            } else if(seqs[1][i] == '-') {
                // deletion;
                weight += no_gap + gap_open;
                state = 1;
            } else {
                // match/mismatch;
                weight +=
                    2 * no_gap + p_marg(seq_pair[0][i - ngap], seq_pair[1][i]);
            }
            break;
        case 1:  // deletion
            if(seqs[0][i] == '-') {
                throw std::runtime_error(
                    "Insertion after deletion is not modeled.");
            } else if(seqs[1][i] == '-') {
                // deletion_ext
                weight += gap_extend;
            } else {
                // match/mismatch
                weight +=
                    gap_stop + p_marg(seq_pair[0][i - ngap], seq_pair[1][i]);
                state = 0;
            }
            break;
        case 2:  // insertion
            if(seqs[0][i] == '-') {
                // insertion_ext
                weight += gap_extend;
                ngap++;
            } else if(seqs[1][i] == '-') {
                // deletion
                weight += gap_stop + gap_open;
                state = 1;
            } else {
                // match/mismatch
                weight += gap_stop + no_gap +
                          p_marg(seq_pair[0][i - ngap], seq_pair[1][i]);
                state = 0;
            }
            break;
        }
    }
    // terminal state weight
    if(state == 0) {
        weight += no_gap;
    } else if(state == 2) {
        weight += gap_stop;
    }
    return weight;
}

/// @private
// GCOVR_EXCL_START
TEST_CASE("alignment_score") {
    coati::alignment_t aln;
    coati::Matrixf P(mg94_p(0.0133, 0.2, {0.308, 0.185, 0.199, 0.308}));
    coati::Matrixf p_marg = marginal_p(P, aln.pi, AmbiguousNucs::AVG);

    aln.data.seqs = {"CTCTGGATAGTG", "CT----ATAGTG"};
    CHECK_EQ(alignment_score(aln, p_marg), doctest::Approx(1.51294f));

    aln.data.seqs = {"CTCT--AT", "CTCTGGAT"};
    CHECK_EQ(alignment_score(aln, p_marg), doctest::Approx(-0.835939f));

    aln.data.seqs = {"ACTCT-A", "ACTCTG-"};
    CHECK_EQ(alignment_score(aln, p_marg), doctest::Approx(-8.73357f));

    aln.data.seqs = {"ACTCTA-", "ACTCTAG"};
    CHECK_EQ(alignment_score(aln, p_marg), doctest::Approx(-0.658564f));
    // different length
    aln.data.seqs = {"CTC", "CT"};
    REQUIRE_THROWS_AS(alignment_score(aln, p_marg), std::invalid_argument);
    // insertion after deletion is not modeled
    aln.data.seqs = {"ATAC-GGGTC", "ATA-GGGGTC"};
    REQUIRE_THROWS_AS(alignment_score(aln, p_marg), std::runtime_error);
}
// GCOVR_EXCL_STOP

/**
 * \brief Sample from a marginal alignment
 *
 * @param[in,out] aln coati::alignment_t alignment data.
 * @param[in] sample_size size_t number of alignments to sample.
 * @param[in] rand coati::random_t random seed generator object.
 *
 */
void marg_sample(coati::alignment_t& aln, size_t sample_size, random_t& rand) {
    coati::Matrixf P(64, 64), p_marg;

    std::ostream* pout(nullptr);
    std::ofstream outfile;
    if(aln.data.out_file.path.empty() || aln.data.out_file.path == "-") {
        pout = &std::cout;
    } else {
        outfile.open(aln.data.out_file.path);
        if(!outfile) {
            throw std::invalid_argument("Opening output file" +
                                        aln.data.out_file.path + "  failed.");
        }
        pout = &outfile;
    }
    std::ostream& out = *pout;

    // check that length of ref sequence is multiple of 3 and gap unit size
    size_t len_a = aln.seq(0).length();
    if((len_a % 3 != 0) || (len_a % aln.gap.len != 0)) {
        throw std::invalid_argument(
            "Length of reference sequence must be multiple of 3.");
    }
    // check that length of descendant sequence is multiple of gap unit size
    if(aln.seq(1).length() % aln.gap.len != 0) {
        throw std::invalid_argument(
            "Length of descendant sequence must be multiple of " +
            std::to_string(aln.gap.len) + ".");
    }

    // encode sequences
    auto anc = aln.seq(0);
    auto des = aln.seq(1);
    coati::utils::sequence_pair_t seq_pair =
        coati::utils::marginal_seq_encoding(anc, des);

    // set substitution matrix according to model
    coati::utils::set_subst(aln);

    // dynamic programming pairwise alignment and traceback
    coati::align_pair_work_t work;
    coati::viterbi(work, seq_pair[0], seq_pair[1], aln);

    out << "[" << std::endl;

    for(size_t i = 0; i < sample_size; ++i) {
        coati::sampleback(work, anc, des, aln, aln.gap.len, rand);

        out << "  {\n    \"aln\": {\n";
        out << "      \"" << aln.name(0) << "\": ";
        out << "\"" << aln.seq(0) << "\",\n";
        out << "      \"" << aln.name(1) << "\": ";
        out << "\"" << aln.seq(1) << "\"\n";
        out << "    },\n";
        out << "    \"weight\": " << ::expf(aln.data.weight) << ",\n";
        out << "    \"log_weight\": " << aln.data.weight << "\n";
        if(i < sample_size - 1) {
            out << "  },";
        } else {
            out << "  }";
        }
        out << std::endl;
    }

    out << "]" << std::endl;
}

/// @private
// GCOVR_EXCL_START
TEST_CASE("marg_sample") {
    // test helper function
    // NOLINTNEXTLINE(misc-unused-parameters)
    auto check_line_eq = [](std::ifstream& in, const std::string_view line) {
        std::string s;
        std::getline(in, s);
        CHECK_EQ(s, line);
    };

    auto test = [check_line_eq](const std::string& seq1,
                                const std::string& seq2,
                                const std::vector<std::string>& expected_s1,
                                const std::vector<std::string>& expected_s2,
                                const std::vector<std::string>& weight,
                                const std::vector<std::string>& lweight) {
        // set seed
        coati::random_t rand;
        const std::vector<std::string> s42{{"42"}};
        auto seed = fragmites::random::string_seed_seq(s42.begin(), s42.end());
        rand.Seed(seed);

        coati::alignment_t aln;
        aln.data = coati::data_t("", {"A", "B"}, {seq1, seq2});
        aln.data.out_file = {{"test-marg_sample.json"}, {".json"}};

        size_t reps{expected_s1.size()};
        coati::marg_sample(aln, reps, rand);

        std::ifstream infile(aln.data.out_file.path);
        REQUIRE(infile.good());

        check_line_eq(infile, "[");
        for(size_t i = 0; i < reps; ++i) {
            check_line_eq(infile, "  {");
            check_line_eq(infile, "    \"aln\": {");
            check_line_eq(infile, "      \"A\": \"" + expected_s1[i] + "\",");
            check_line_eq(infile, "      \"B\": \"" + expected_s2[i] + "\"");
            check_line_eq(infile, "    },");
            check_line_eq(infile, "    \"weight\": " + weight[i] + ",");
            check_line_eq(infile, "    \"log_weight\": " + lweight[i]);
            if(i < reps - 1) {
                check_line_eq(infile, "  },");
            } else {
                check_line_eq(infile, "  }");
            }
        }
        check_line_eq(infile, "]");
        infile.close();
        REQUIRE(std::filesystem::remove(aln.data.out_file.path));
    };

    SUBCASE("sample size 1") {
        std::vector<std::string> seq1{"CC--CCCC"};
        std::vector<std::string> seq2{"CCCCCCCC"};
        std::vector<std::string> weight{"0.031239"};
        std::vector<std::string> lweight{"-3.46609"};
        test("CCCCCC", "CCCCCCCC", seq1, seq2, weight, lweight);
    }
    SUBCASE("sample size 3") {
        std::vector<std::string> seq1{"CC--CCCC", "CCCCCC--", "CCCCC--C"};
        std::vector<std::string> seq2{"CCCCCCCC", "CCCCCCCC", "CCCCCCCC"};
        std::vector<std::string> weight{"0.031239", "0.499854", "0.249923"};
        std::vector<std::string> lweight{"-3.46609", "-0.69344", "-1.3866"};
        test("CCCCCC", "CCCCCCCC", seq1, seq2, weight, lweight);
    }
    SUBCASE("length of reference not multiple of 3") {
        coati::random_t rand;
        coati::alignment_t aln;
        aln.data = coati::data_t("", {"A", "B"}, {"C", "CCC"});
        CHECK_THROWS_AS(marg_sample(aln, 1, rand), std::invalid_argument);
    }
    SUBCASE("length of descendant no multiple of gap len") {
        coati::random_t rand;
        coati::alignment_t aln;
        aln.data = coati::data_t("", {"A", "B"}, {"CCC", "CCCC"});
        aln.gap.len = 3;
        std::vector<std::string> seq1, seq2, weight, lweight;
        CHECK_THROWS_AS(marg_sample(aln, 1, rand), std::invalid_argument);
    }
    SUBCASE("error opening output file") {
        coati::random_t rand;
        coati::alignment_t aln;
        aln.data.out_file = {{"."}, {".fasta"}};
        CHECK_THROWS_AS(marg_sample(aln, 1, rand), std::invalid_argument);
    }
}
// GCOVR_EXCL_STOP

}  // namespace coati
