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
#include <coati/io.hpp>
#include <coati/json.hpp>
#include <coati/mg94q.tcc>
#include <coati/mutation_coati.hpp>
#include <filesystem>

namespace coati {

/**
 * @brief Pairwise alignment using dynamic programming and a marginal model.
 *
 * Alignment of two sequences via dynamic programming using an affine
 *  (geometrical) gap model and a marginal codon substitution model.
 *
 * @param[in,out] aln coati::alignment_t alignment information.
 *
 * @retval true successful run.
 */
bool marg_alignment(coati::alignment_t& aln) {
    std::ofstream out_w;

    // read input data
    aln.data = coati::io::read_input(aln);

    // set substitution matrix according to model
    coati::utils::set_subst(aln);

    // if -s or --score, score alignment and exit
    if(aln.score) {
        std::cout << alignment_score(aln, aln.subst_matrix) << std::endl;
        return true;
    }

    // process input sequences
    coati::utils::process_marginal(aln);

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
        std::cerr << "ERROR: sequences to align exceed available memory."
                  << std::endl;
        return EXIT_FAILURE;
    } catch(const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << std::endl;
        return false;
    }
    coati::traceback_viterbi(work, anc, des, aln, aln.gap.len);

    // handle end stop codons
    coati::utils::restore_end_stops(aln.data, aln.gap);

    // write alignment
    coati::io::write_output(aln);
    return true;
}

/// @private
// GCOVR_EXCL_START
TEST_CASE("marg_alignment") {
    // NOLINTNEXTLINE(misc-unused-parameters)
    auto test_fasta = [](alignment_t& aln, const data_t& expected,
                         const std::string& file) {
        std::ofstream out;
        out.open(aln.data.path);
        REQUIRE(out);
        out << file;
        out.close();

        REQUIRE(marg_alignment(aln));

        std::ifstream infile(aln.output);
        std::string s1, s2;

        infile >> s1 >> s2;
        CHECK_EQ(s1, expected.names[0]);
        CHECK_EQ(s2, expected.seqs[0]);

        infile >> s1 >> s2;
        CHECK_EQ(s1, expected.names[1]);
        CHECK_EQ(s2, expected.seqs[1]);
        REQUIRE(std::filesystem::remove(aln.output));
        REQUIRE(std::filesystem::remove("test-marg.fasta"));
    };

    // NOLINTNEXTLINE(misc-unused-parameters)
    auto test_phylip = [](alignment_t& aln, const data_t& expected,
                          const std::string& file) {
        std::ofstream out;
        out.open(aln.data.path);
        REQUIRE(out);
        out << file;
        out.close();

        REQUIRE(marg_alignment(aln));

        std::ifstream infile(aln.output);
        std::string s1, s2;

        infile >> s1 >> s2;
        CHECK_EQ(s1, "2");
        CHECK_EQ(std::stoul(s2), aln.seq(0).size());

        infile >> s1 >> s2;
        CHECK_EQ(s1, expected.names[0]);
        CHECK_EQ(s2, expected.seqs[0]);

        infile >> s1 >> s2;
        CHECK_EQ(s1, expected.names[1]);
        CHECK_EQ(s2, expected.seqs[1]);
        REQUIRE(std::filesystem::remove(aln.output));
        REQUIRE(std::filesystem::remove("test-marg.fasta"));
    };

    alignment_t aln;
    data_t expected;
    SUBCASE("Alignment - output fasta") {
        std::string file{">1\nCTCTGGATAGTG\n>2\nCTATAGTG\n"};
        aln.data.path = "test-marg.fasta";
        aln.model = "mar-mg";
        aln.output = "test-marg_alignment-fasta.fasta";
        expected = data_t("", {">1", ">2"}, {"CTCTGGATAGTG", "CT----ATAGTG"});

        test_fasta(aln, expected, file);
    }
    SUBCASE("Alignment - output fasta - use aln.refs") {
        std::string file{">1\nCTATAGTG\n>2\nCTCTGGATAGTG\n"};
        aln.data.path = "test-marg.fasta";
        aln.model = "mar-mg";
        aln.output = "test-marg_alignment-fasta.fasta";
        aln.refs = "2";

        expected = data_t("", {">2", ">1"}, {"CTCTGGATAGTG", "CT----ATAGTG"});

        test_fasta(aln, expected, file);
    }
    SUBCASE("Alignment - output fasta - use aln.refs - no change") {
        std::string file{">1\nCTCTGGATAGTG\n>2\nCTATAGTG\n"};
        aln.data.path = "test-marg.fasta";
        aln.model = "mar-mg";
        aln.output = "test-marg_alignment-fasta.fasta";
        aln.refs = "1";

        expected = data_t("", {">1", ">2"}, {"CTCTGGATAGTG", "CT----ATAGTG"});

        test_fasta(aln, expected, file);
    }
    SUBCASE("Alignment - output phylip") {
        std::string file{">1\nGCGACTGTT\n>2\nGCGATTGCTGTT\n"};
        aln.data.path = "test-marg.fasta";
        aln.model = "mar-mg";
        aln.output = "test-marg_alignment-phylip.phy";

        expected = data_t("", {"1", "2"}, {"GCGA---CTGTT", "GCGATTGCTGTT"});

        test_phylip(aln, expected, file);
    }
    SUBCASE("Alignment - output phylip - use aln.refn") {
        std::string file{">A\nGCGATTGCTGTT\n>B\nGCGACTGTT\n"};
        aln.data.path = "test-marg.fasta";
        aln.model = "mar-ecm";
        aln.output = "test-marg_alignment-phylip.phy";
        aln.rev = true;

        expected = data_t("", {"B", "A"}, {"GCGA---CTGTT", "GCGATTGCTGTT"});

        test_phylip(aln, expected, file);
    }
    SUBCASE("Alignment 2 dels - output phylip") {
        std::string file{">1\nACGTTAAGGGGT\n>2\nACGAAT\n"};
        aln.data.path = "test-marg.fasta";
        aln.model = "mar-mg";
        aln.output = "test-marg_alignment-phylip2.phy";

        expected = data_t("", {"1", "2"}, {"ACGTTAAGGGGT", "ACG--AA----T"});

        test_phylip(aln, expected, file);
    }
    SUBCASE("Alignment with gap length multiple of 3") {
        std::string file{">1\nACGTTAAGGGGT\n>2\nACGAAT\n"};
        aln.data.path = "test-marg.fasta";
        aln.model = "mar-mg";
        aln.output = "test-marg_alignment-no-frameshifts.fa";
        aln.gap.len = 3;

        expected = data_t("", {">1", ">2"}, {"ACGTTAAGGGGT", "AC------GAAT"});

        test_fasta(aln, expected, file);
    }
    SUBCASE("Alignment with ambiguous nucleotides - SUM") {
        std::string file{">1\nCTCTGGATAGTG\n>2\nCTATAGTR\n"};
        aln.data.path = "test-marg.fasta";
        aln.model = "mar-mg";
        aln.output = "test-marg_alignment_ambiguous.fa";

        expected = data_t("", {">1", ">2"}, {"CTCTGGATAGTG", "CT----ATAGTR"});
        test_fasta(aln, expected, file);
    }
    SUBCASE("Alignment with ambiguous nucleotides - BEST") {
        std::string file{">1\nCTCTGGATAGTG\n>2\nCTATAGTR\n"};
        aln.data.path = "test-marg.fasta";
        aln.model = "mar-mg";
        aln.output = "test-marg_alignment_ambiguous.fa";
        aln.amb = coati::AmbiguousNucs::BEST;

        expected = data_t("", {">1", ">2"}, {"CTCTGGATAGTG", "CT----ATAGTR"});
        test_fasta(aln, expected, file);
    }
    SUBCASE("Alignment with gap length multiple of 3 - fail") {
        std::ofstream out;
        out.open("test-marg.fasta");
        out << ">1\nGCGATTGCTGT\n>2\nGCGACTGTT\n";
        out.close();
        aln.model = "mar-mg";
        aln.data.path = "test-marg.fasta";
        aln.output = "test-marg_alignment-no-frameshifts-f.fasta";
        aln.gap.len = 3;

        REQUIRE_THROWS_AS(marg_alignment(aln), std::invalid_argument);
        CHECK(std::filesystem::remove("test-marg.fasta"));
    }
    SUBCASE("Length of descendant not multiple of gap.len") {
        std::ofstream out;
        out.open("test-marg.fasta");
        out << ">A\nCTCGGA\n>B\nCTCGG\n";
        out.close();
        aln.data.path = "test-marg.fasta";
        aln.model = "mar-mg";
        aln.gap.len = 3;
        REQUIRE_THROWS_AS(marg_alignment(aln), std::invalid_argument);
        CHECK(std::filesystem::remove("test-marg.fasta"));
    }
    SUBCASE("Score alignment") {
        std::ofstream out;
        out.open("test-marg.fasta");
        out << ">1\nCTCTGGATAGTG\n>2\nCT----ATAGTG\n";
        out.close();
        aln.data.path = "test-marg.fasta";
        aln.model = "mar-mg";
        aln.output = "test-marg_alignment-score.fasta";
        aln.score = true;

        REQUIRE(marg_alignment(aln));
        CHECK(std::filesystem::remove("test-marg.fasta"));
    }
    SUBCASE("Score alignment diff length - fail") {
        std::ofstream out;
        out.open("test-marg.fasta");
        out << ">1\nCTCTGGATAGTG\n>2\nCTATAGTG\n";
        out.close();
        aln.data.path = "test-marg.fasta";
        aln.model = "mar-mg";
        aln.output = "test-marg_alignment-score-f.fasta";
        aln.score = true;

        coati::data_t result(aln.output);

        REQUIRE_THROWS_AS(marg_alignment(aln), std::invalid_argument);
        CHECK(std::filesystem::remove("test-marg.fasta"));
    }
    SUBCASE("Name of reference not found") {
        std::ofstream out;
        out.open("test-marg.fasta");
        out << ">1\nCTCTGGATAGTG\n>2\nCTATAGTG\n";
        out.close();
        aln.data.path = "test-marg.fasta";
        aln.refs = "seq_name";

        REQUIRE_THROWS_AS(marg_alignment(aln), std::invalid_argument);
        CHECK(std::filesystem::remove("test-marg.fasta"));
    }
    SUBCASE("User-provided codon substitution matrix") {
        std::ofstream outfile;
        coati::Matrix<coati::float_t> P(
            mg94_p(0.0133, 0.2, {0.308, 0.185, 0.199, 0.308}));

        const std::vector<std::string> codons = {
            "AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA",
            "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT", "CAA", "CAC",
            "CAG", "CAT", "CCA", "CCC", "CCG", "CCT", "CGA", "CGC", "CGG",
            "CGT", "CTA", "CTC", "CTG", "CTT", "GAA", "GAC", "GAG", "GAT",
            "GCA", "GCC", "GCG", "GCT", "GGA", "GGC", "GGG", "GGT", "GTA",
            "GTC", "GTG", "GTT", "TAC", "TAT", "TCA", "TCC", "TCG", "TCT",
            "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT"};

        outfile.open("test-marg-matrix.csv");
        REQUIRE(outfile);

        float Q[3721]{0.0f};
        for(auto i = 0; i < 587; i++) {
            Q[mg94_indexes[i]] = mg94Q[i];
        }

        outfile << "0.0133" << std::endl;  // branch length
        for(auto i = 0; i < 61; i++) {
            for(auto j = 0; j < 61; j++) {
                outfile << codons[i] << "," << codons[j] << "," << Q[i * 61 + j]
                        << std::endl;
            }
        }

        outfile.close();
        std::string file{">1\nCTCTGGATAGTG\n>2\nCTATAGTG\n"};
        aln.data.path = "test-marg.fasta";
        aln.rate = "test-marg-matrix.csv";
        aln.output = "test-marg_alignment-fasta.fasta";

        expected = data_t("", {">1", ">2"}, {"CTCTGGATAGTG", "CT----ATAGTG"});

        test_fasta(aln, expected, file);
        REQUIRE(std::filesystem::remove("test-marg-matrix.csv"));
    }
    SUBCASE("Numb of sequences != 2 - fail") {
        aln.data.path = "test-marg.fasta";
        std::ofstream out;
        out.open(aln.data.path);
        REQUIRE(out);
        out << ">1\nCTCTGGATAGTG\n";
        out.close();
        CHECK_THROWS_AS(marg_alignment(aln), std::invalid_argument);

        aln.data.path = "test-marg.fasta";
        out.open("test-marg.fasta");
        REQUIRE(out);
        out << ">1\nCTCTGGATAGTG\n>2\nCTATAGTG\n>3\nCTCTGGGTG\n";
        out.close();
        CHECK_THROWS_AS(marg_alignment(aln), std::invalid_argument);
        REQUIRE(std::filesystem::remove("test-marg.fasta"));
    }
}
// GCOVR_EXCL_STOP

/**
 * @brief Score alignment using marginal model.
 *
 * @param[in] aln coati::alignment_t input parameters.
 * @param[in] p_marg coati::Matrixf substitution matrix.
 *
 * @retval float alignment score.
 */
float alignment_score(coati::alignment_t& aln, const coati::Matrixf& p_marg) {
    // Remove gaps from alignment and return an expanded cigar string of
    // alignment
    std::string cigar = coati::utils::process_alignment(aln);

    coati::utils::sequence_pair_t seq_pair =
        coati::utils::marginal_seq_encoding(aln.data.seqs[0], aln.data.seqs[1]);

    coati::semiring::tropical trig;  // tropical rig
    // calculate log(1-g) log(1-e) log(g) log(e) log(pi)
    float_t no_gap = trig.from_linear_1mf(aln.gap.open);
    float_t gap_stop = trig.from_linear_1mf(aln.gap.extend);
    float_t gap_open = trig.from_linearf(aln.gap.open);
    float_t gap_extend = trig.from_linearf(aln.gap.extend);
    std::vector<coati::float_t> pi{aln.pi};
    std::transform(pi.cbegin(), pi.cend(), pi.begin(),
                   [](auto value) { return ::logf(value); });

    enum struct State { MATCH, GAP };
    auto state = State::MATCH;
    float score{0.f};
    size_t nins{0}, ndel{0};
    size_t apos = 0, bpos = 0;

    for(size_t i = 0; i < cigar.length(); i++) {
        switch(state) {
        case State::MATCH:
            if(cigar[i] == 'I') {  // insertion;
                nins++;
                bpos++;
                state = State::GAP;
            } else if(cigar[i] == 'D') {  // deletion;
                ndel++;
                apos++;
                state = State::GAP;
            } else {  // match/mismatch;
                score =
                    trig.times(score, no_gap, no_gap,
                               p_marg(seq_pair[0][apos], seq_pair[1][bpos]));
                apos++;
                bpos++;
            }
            break;
        case State::GAP:
            if(cigar[i] == 'I') {  // insertion_extension
                nins++;
                bpos++;
            } else if(cigar[i] == 'D') {  // deletion_ext
                ndel++;
                apos++;
            } else {  // match/mismatch
                assert(nins > 0 || ndel > 0);
                if(nins == 0) {  // score deletions
                    score =
                        trig.times(score, no_gap, gap_open,
                                   trig.power(gap_extend, ndel - 1), gap_stop);
                } else if(ndel == 0) {  // score insertions
                    score = trig.times(score, gap_open,
                                       trig.power(gap_extend, nins - 1),
                                       gap_stop, no_gap);
                } else {  // score both insertions and deletions
                    score = trig.times(score, gap_open, gap_open,
                                       trig.power(gap_extend, nins + ndel - 2),
                                       gap_stop, gap_stop);
                }
                score = trig.times(
                    score, p_marg(seq_pair[0][apos], seq_pair[1][bpos]));
                nins = ndel = 0;
                state = State::MATCH;
                apos++;
                bpos++;
            }
            break;
        }
    }
    assert(apos == seq_pair[0].length());
    assert(bpos == seq_pair[1].length());
    // terminal state score
    if(state == State::MATCH) {
        score = trig.times(score, no_gap, no_gap);
    } else if(state == State::GAP) {
        if(nins == 0) {  // score deletions
            score = trig.times(score, no_gap, gap_open,
                               trig.power(gap_extend, ndel - 1), gap_stop);
        } else if(ndel == 0) {  // score insertions
            score =
                trig.times(score, gap_open, trig.power(gap_extend, nins - 1),
                           gap_stop, no_gap);
        } else {  // score both insertions and deletions
            score = trig.times(score, gap_open, gap_open,
                               trig.power(gap_extend, nins + ndel - 2),
                               gap_stop, gap_stop, no_gap);
        }
    }

    // handle end stop codons
    aln.data.score = score;
    coati::utils::restore_end_stops(aln.data, aln.gap);

    return aln.data.score;
}

/// @private
// GCOVR_EXCL_START
TEST_CASE("alignment_score") {
    // NOLINTNEXTLINE(misc-unused-parameters)
    auto test = [](const std::string& anc, const std::string& des, float exp) {
        coati::alignment_t aln;
        coati::Matrixf P(mg94_p(0.0133, 0.2, {0.308, 0.185, 0.199, 0.308}));
        coati::Matrixf p_marg =
            marginal_p(P, aln.pi, AmbiguousNucs::SUM, MarginalSubst::SUM);
        aln.data.names = {"A", "B"};
        aln.data.seqs = {anc, des};
        CHECK_EQ(alignment_score(aln, p_marg), doctest::Approx(exp));
    };

    test("CTCTGGATAGTG", "CT----ATAGTG", 1.50914f);
    test("CTCT--AT", "CTCTGGAT", -0.83906f);
    test("ACTCT-A", "ACTCTG-", -10.52864);
    test("ATGCTTTAC", "ATGCT-TAC", 2.13593f);
    test("ATGCTT---", "ATGCTTTGA", 0.70607f);
    test("A-CTAAC", "ACCTAAG", -8.2786f);
    test("ACT---", "ACTCTG", -5.04197);
    test("ACTCTA", "ACT---", -5.04197);
    test("ACT----", "ACT-CTG", -5.04197);
    test("AAAAAA---AAA", "AAA---AAAAAA", -11.09557);
    test("AAA---AAAAAA", "AAAAAA---AAA", -11.09557);
    test("AAA-A-A-AAAA", "AAAA-A-A-AAA", -11.09557);
    test("---AAAAAA", "AAAAAAAAA", -2.03242);
    test("AAAAAA---", "AAAAAAAAA", -2.03242);
    test("AAAAAAAAA", "---AAAAAA", -2.03242);
    test("AAAAAAAAA", "AAAAAA---", -2.03242);
    test("ACTCTA", "ACTC--", -3.18537f);
    // Theses alignment will be scored as ACTCTA--- / ACTC--TAG
    test("ACTCTA-", "ACTCTAG", -10.45777f);
    test("ACTCTA--", "ACTCT-AG", -10.45777f);

    // NOLINTNEXTLINE(misc-unused-parameters)
    auto test_fail = [](const std::string& anc, const std::string& des) {
        coati::alignment_t aln;
        coati::Matrixf P(mg94_p(0.0133, 0.2, {0.308, 0.185, 0.199, 0.308}));
        coati::Matrixf p_marg =
            marginal_p(P, aln.pi, AmbiguousNucs::SUM, MarginalSubst::SUM);
        aln.data.names = {"A", "B"};
        aln.data.seqs = {anc, des};
        CHECK_THROWS_AS(alignment_score(aln, p_marg), std::invalid_argument);
    };

    // more than 2 sequences
    test_fail("ATACGGGTC", "");
    // length of ref is not multiple of 3
    test_fail("ATAC", "ATA-");
}
// GCOVR_EXCL_STOP

/**
 * @brief Sample from a marginal alignment.
 *
 * @param[in,out] aln coati::alignment_t alignment data.
 * @param[in] sample_size size_t number of alignments to sample.
 * @param[in] rand coati::random_t random seed generator object.
 *
 */
void marg_sample(coati::alignment_t& aln, size_t sample_size, random_t& rand) {
    coati::Matrixf P(61, 61), p_marg;

    // read input data
    aln.data = coati::io::read_input(aln);
    if(aln.data.size() != 2) {
        throw std::invalid_argument("Exactly two sequences required.");
    }

    // set output pointer
    std::ostream* pout(nullptr);
    std::ofstream outfile;
    if(aln.output.empty() || aln.output == "-") {
        pout = &std::cout;
    } else {
        outfile.open(aln.output);
        if(!outfile) {
            throw std::invalid_argument("Opening output file " +
                                        aln.output.string() + " failed.");
        }
        pout = &outfile;
    }
    std::ostream& out = *pout;

    // check that length of ref sequence is multiple of 3 and gap unit size
    size_t len_a = aln.seq(0).length();
    if(len_a % 3 != 0 || len_a % aln.gap.len != 0) {
        throw std::invalid_argument(
            "Length of reference sequence must be multiple of 3.");
    }
    // check that length of descendant sequence is multiple of gap unit size
    if(aln.seq(1).length() % aln.gap.len != 0) {
        throw std::invalid_argument(
            "Length of descendant sequence must be multiple of " +
            std::to_string(aln.gap.len) + ".");
    }

    coati::utils::trim_end_stops(aln.data);

    // encode sequences
    auto anc = aln.seq(0);
    auto des = aln.seq(1);
    coati::utils::sequence_pair_t seq_pair =
        coati::utils::marginal_seq_encoding(anc, des);

    // set substitution matrix according to model
    coati::utils::set_subst(aln);

    // dynamic programming pairwise alignment and sampleback
    coati::align_pair_work_t work;
    coati::forward(work, seq_pair[0], seq_pair[1], aln);

    // sample and print as many aligments as required (sample_size)
    for(size_t i = 0; i < sample_size; ++i) {
        coati::sampleback(work, anc, des, aln, aln.gap.len, rand);
        coati::utils::restore_end_stops(aln.data, aln.gap);
        write_json(aln.data, out, i, sample_size);
    }
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
                                const std::vector<std::string>& lscore) {
        // set seed
        coati::random_t rand;
        const std::vector<std::string> s42{{"42"}};
        auto seed = fragmites::random::string_seed_seq(s42.begin(), s42.end());
        rand.Seed(seed);

        coati::alignment_t aln;
        aln.data.path = "test-marg.fasta";
        aln.output = "test-marg_sample.json";
        std::ofstream out;
        out.open(aln.data.path);
        REQUIRE(out);
        out << ">A\n" << seq1 << "\n>B\n" << seq2 << std::endl;
        out.close();

        size_t reps{expected_s1.size()};
        coati::marg_sample(aln, reps, rand);

        std::ifstream infile(aln.output);
        REQUIRE(infile.good());

        check_line_eq(infile, "[");
        for(size_t i = 0; i < reps; ++i) {
            check_line_eq(infile, "{");
            check_line_eq(infile, "  \"alignment\": {");
            check_line_eq(infile, R"(    "A": ")" + expected_s1[i] + "\",");
            check_line_eq(infile, R"(    "B": ")" + expected_s2[i] + "\"");
            check_line_eq(infile, "  },");
            check_line_eq(infile, "  \"score\": " + lscore[i]);
            if(i < reps - 1) {
                check_line_eq(infile, "},");
            } else {
                check_line_eq(infile, "}");
            }
        }
        check_line_eq(infile, "]");
        infile.close();
        REQUIRE(std::filesystem::remove(aln.output));
        REQUIRE(std::filesystem::remove("test-marg.fasta"));
    };

    SUBCASE("sample size 1") {
        std::vector<std::string> seq1{"CC--CCCC"};
        std::vector<std::string> seq2{"CCCCCCCC"};
        std::vector<std::string> lscore{"-1.9466571807861328"};
        test("CCCCCC", "CCCCCCCC", seq1, seq2, lscore);
    }
    SUBCASE("sample size 1 - deletion") {
        std::vector<std::string> seq1{"CCCCCC"};
        std::vector<std::string> seq2{"--CCCC"};
        std::vector<std::string> lscore{"-1.6172490119934082"};
        test("CCCCCC", "CCCC", seq1, seq2, lscore);
    }
    SUBCASE("sample size 3") {
        std::vector<std::string> seq1{"CC--CCCC", "CCCCCC--", "CCCC--CC"};
        std::vector<std::string> seq2{"CCCCCCCC", "CCCCCCCC", "CCCCCCCC"};
        std::vector<std::string> lscore{"-1.9466571807861328",
                                        "-1.9466569423675537",
                                        "-1.9466572999954224"};
        test("CCCCCC", "CCCCCCCC", seq1, seq2, lscore);
    }
    SUBCASE("length of reference not multiple of 3") {
        coati::random_t rand;
        coati::alignment_t aln;
        aln.data.path = "marg-sample.fasta";
        std::ofstream outfile;
        outfile.open(aln.data.path);
        REQUIRE(outfile);
        outfile << ">seq1\nAC\n>seq2\nACG\n";
        outfile.close();
        CHECK_THROWS_AS(marg_sample(aln, 1, rand), std::invalid_argument);
        REQUIRE(std::filesystem::remove("marg-sample.fasta"));
    }
    SUBCASE("length of descendant no multiple of gap len") {
        coati::random_t rand;
        coati::alignment_t aln;
        aln.data.path = "marg-sample.fasta";
        aln.gap.len = 3;
        std::ofstream outfile;
        outfile.open(aln.data.path);
        REQUIRE(outfile);
        outfile << ">A\nCCC\n>B\nCCCC\n";
        outfile.close();
        CHECK_THROWS_AS(marg_sample(aln, 1, rand), std::invalid_argument);
        REQUIRE(std::filesystem::remove("marg-sample.fasta"));
    }
    SUBCASE("error opening output file") {
        coati::random_t rand;
        coati::alignment_t aln;
        aln.data.path = "marg-sample.fasta";
        aln.output = "..fasta";
        CHECK_THROWS_AS(marg_sample(aln, 1, rand), std::invalid_argument);
    }
    SUBCASE("Number of seqs  != 2") {
        coati::random_t rand;
        coati::alignment_t aln;
        aln.data.path = "marg-sample.fasta";
        std::ofstream outfile;
        outfile.open(aln.data.path);
        REQUIRE(outfile);
        outfile << ">A\nCCC\n";
        outfile.close();
        CHECK_THROWS_AS(marg_sample(aln, 1, rand), std::invalid_argument);
        aln.data.path = "marg-sample.fasta";
        outfile.open("marg-sample.fasta");
        REQUIRE(outfile);
        outfile << ">A\nCCC\n>B\nCCC\n>C\nCCC\n";
        outfile.close();
        CHECK_THROWS_AS(marg_sample(aln, 1, rand), std::invalid_argument);
        REQUIRE(std::filesystem::remove("marg-sample.fasta"));
    }
}
// GCOVR_EXCL_STOP

}  // namespace coati
