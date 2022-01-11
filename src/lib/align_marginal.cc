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

    // score alignment
    if(aln.score) {
        std::cout << alignment_score(aln, aln.subst_matrix) << std::endl;
        return true;
    }

    // check that length of ref sequence is multiple of 3 and gap unit size
    size_t len_a = aln.data.seqs[0].length();
    if((len_a % 3 != 0) && (len_a % aln.gap.len != 0)) {
        throw std::invalid_argument(
            "Length of reference sequence must be multiple of 3.");
    }
    // check that length of descendant sequence is multiple of gap unit size
    if(aln.data.seqs[1].length() % aln.gap.len != 0) {
        throw std::invalid_argument(
            "Length of descendant sequence must be multiple of " +
            std::to_string(aln.gap.len) + ".");
    }

    // encode sequences
    auto anc = aln.data.seqs[0];
    auto des = aln.data.seqs[1];
    coati::utils::sequence_pair_t seq_pair =
        coati::utils::marginal_seq_encoding(anc, des);

    // dynamic programming pairwise alignment and traceback
    coati::align_pair_work_mem_t work;
    coati::viterbi_mem(work, seq_pair[0], seq_pair[1], aln);
    coati::traceback(work, anc, des, aln, aln.gap.len);

    if(!aln.weight_file.empty()) {  // save weight and filename
        out_w.open(aln.weight_file, std::ios::app | std::ios::out);
        out_w << aln.data.path << "," << aln.model << "," << aln.data.weight
              << std::endl;
        out_w.close();
    }

    // write alignment
    return coati::utils::write_output(aln.data);
}

/// @private
TEST_CASE("marg_alignment") {
    SUBCASE("Alignment - output fasta") {
        coati::alignment_t aln;
        aln.data = coati::data_t("", {"1", "2"}, {"CTCTGGATAGTG", "CTATAGTG"});
        aln.model = "m-coati";
        aln.weight_file = "score.log";
        aln.data.out_file = {{"test-marg_alignment-fasta.fasta"}, {".fasta"}};

        coati::utils::set_subst(aln);

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
        CHECK(s1 == ">1");
        CHECK(s2 == "CTCTGGATAGTG");

        infile >> s1 >> s2;
        CHECK(s1 == ">2");
        CHECK(s2 == "CT----ATAGTG");
        CHECK(std::filesystem::remove(aln.data.out_file.path));

        std::ifstream inweight(aln.weight_file);
        std::string s;
        inweight >> s;
        CHECK(std::filesystem::remove(aln.weight_file));
        CHECK(s.substr(s.length() - 7) == "1.51294");
    }

    SUBCASE("Alignment - output phylip") {
        coati::alignment_t aln;
        aln.data = coati::data_t("", {"1", "2"}, {"GCGACTGTT", "GCGATTGCTGTT"});
        aln.model = "m-coati";
        aln.weight_file = "";
        aln.data.out_file = {{"test-marg_alignment-phylip.phy"}, {".phy"}};

        coati::utils::set_subst(aln);

        if(std::filesystem::exists(aln.data.out_file.path)) {
            std::filesystem::remove(aln.data.out_file.path);
        }

        REQUIRE(marg_alignment(aln));

        std::ifstream infile(aln.data.out_file.path);
        std::string s1, s2;

        infile >> s1 >> s2;
        CHECK(s1 == "2");
        CHECK(s2 == "12");

        infile >> s1 >> s2;
        CHECK(s1 == "1");
        CHECK(s2 == "GCGA---CTGTT");

        infile >> s1 >> s2;
        CHECK(s1 == "2");
        CHECK(s2 == "GCGATTGCTGTT");

        CHECK(std::filesystem::remove(aln.data.out_file.path));
    }

    SUBCASE("Alignment 2 dels - output phylip") {
        coati::alignment_t aln;
        aln.data = coati::data_t("", {"1", "2"}, {"ACGTTAAGGGGT", "ACGAAT"});
        aln.model = "m-coati";
        aln.weight_file = "";
        aln.data.out_file = {{"test-marg_alignment-phylip2.phy"}, {".phy"}};

        coati::utils::set_subst(aln);

        if(std::filesystem::exists(aln.data.out_file.path)) {
            std::filesystem::remove(aln.data.out_file.path);
        }

        REQUIRE(marg_alignment(aln));

        std::ifstream infile(aln.data.out_file.path);
        std::string s1, s2;

        infile >> s1 >> s2;
        CHECK(s1 == "2");
        CHECK(s2 == "12");

        infile >> s1 >> s2;
        CHECK(s1 == "1");
        CHECK(s2 == "ACGTTAAGGGGT");

        infile >> s1 >> s2;
        CHECK(s1 == "2");
        CHECK(s2 == "ACG--AA----T");

        CHECK(std::filesystem::remove(aln.data.out_file.path));
    }

    SUBCASE("Alignment with gap length multiple of 3") {
        coati::alignment_t aln;
        aln.data = coati::data_t("", {"1", "2"}, {"ACGTTAAGGGGT", "ACGAAT"});
        aln.model = "m-coati";
        aln.weight_file = "";
        aln.data.out_file = {{"test-marg_alignment-no-frameshifts.fa"},
                             {".fa"}};
        aln.gap.len = 3;

        coati::utils::set_subst(aln);

        if(std::filesystem::exists(aln.data.out_file.path)) {
            std::filesystem::remove(aln.data.out_file.path);
        }

        REQUIRE(marg_alignment(aln));
        std::ifstream infile(aln.data.out_file.path);  // input file
        std::string s1, s2;

        infile >> s1 >> s2;
        CHECK(s1 == ">1");
        CHECK(s2 == "ACG---TTAAGGGGT");

        infile >> s1 >> s2;
        CHECK(s1 == ">2");
        CHECK(s2 == "ACGAAT---------");

        CHECK(std::filesystem::remove(aln.data.out_file.path));
    }

    SUBCASE("Alignment with gap length multiple of 3 - fail") {
        coati::alignment_t aln;
        aln.data = coati::data_t("", {"1", "2"}, {"GCGATTGCTGT", "GCGACTGTT"});
        aln.model = "m-coati";
        aln.weight_file = "";
        aln.data.out_file = {{"test-marg_alignment-no-frameshifts-f.fasta"},
                             {".fasta"}};
        aln.gap.len = 3;

        coati::utils::set_subst(aln);
        REQUIRE_THROWS_AS(marg_alignment(aln), std::invalid_argument);
    }

    SUBCASE("Score alignment") {
        coati::alignment_t aln;
        aln.data =
            coati::data_t("", {"1", "2"}, {"CTCTGGATAGTG", "CT----ATAGTG"});
        aln.model = "m-coati";
        aln.weight_file = "";
        aln.data.out_file = {{"test-marg_alignment-score.fasta"}, {".fasta"}};
        aln.score = true;

        coati::utils::set_subst(aln);
        REQUIRE(marg_alignment(aln));
    }

    SUBCASE("Score alignment - fail") {
        coati::alignment_t aln;
        aln.data = coati::data_t("", {"1", "2"}, {"CTCTGGATAGTG", "CTATAGTG"});
        aln.model = "m-coati";
        aln.weight_file = "";
        aln.data.out_file = {{"test-marg_alignment-score-f.fasta"}, {".fasta"}};
        aln.score = true;

        coati::data_t result(aln.data.out_file.path);
        coati::utils::set_subst(aln);

        REQUIRE_THROWS_AS(marg_alignment(aln), std::invalid_argument);
    }
}

/**
 * \brief Score alignment using marginal model.
 *
 * @param[in] args coati::args_t input parameters.
 * @param[in] p_marg coati::Matrixf substitution matrix.
 *
 * \return alignment score (float).
 */
float alignment_score(coati::alignment_t& aln, coati::Matrixf& p_marg) {
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
    for(size_t i = 0; i < 4; i++) {
        pi[i] = ::logf(pi[i]);
    }

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
TEST_CASE("alignment_score") {
    coati::alignment_t aln;
    coati::Matrixf P(mg94_p(0.0133, 0.2, {0.308, 0.185, 0.199, 0.308}));
    coati::Matrixf p_marg = marginal_p(P, aln.pi);

    aln.data.seqs = {"CTCTGGATAGTG", "CT----ATAGTG"};
    REQUIRE(alignment_score(aln, p_marg) == doctest::Approx(1.51294f));

    aln.data.seqs = {"CTCT--AT", "CTCTGGAT"};
    REQUIRE(alignment_score(aln, p_marg) == doctest::Approx(-0.835939f));

    aln.data.seqs = {"CTC", "CT"};
    REQUIRE_THROWS_AS(alignment_score(aln, p_marg), std::invalid_argument);
}

void marg_sample(coati::alignment_t& aln, size_t sample_size, random_t& rand) {
    coati::Matrixf P(64, 64), p_marg;

    std::ostream* pout;
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
    size_t len_a = aln.data.seqs[0].length();
    if((len_a % 3 != 0) && (len_a % aln.gap.len != 0)) {
        throw std::invalid_argument(
            "Length of reference sequence must be multiple of 3.");
    }
    // check that length of descendant sequence is multiple of gap unit size
    if(aln.data.seqs[1].length() % aln.gap.len != 0) {
        throw std::invalid_argument(
            "Length of descendant sequence must be multiple of " +
            std::to_string(aln.gap.len) + ".");
    }

    // encode sequences
    auto anc = aln.data.seqs[0];
    auto des = aln.data.seqs[1];
    coati::utils::sequence_pair_t seq_pair =
        coati::utils::marginal_seq_encoding(anc, des);

    // dynamic programming pairwise alignment and traceback
    coati::align_pair_work_t work;
    coati::viterbi(work, seq_pair[0], seq_pair[1], aln);

    out << "[" << std::endl;

    for(size_t i = 0; i < sample_size; ++i) {
        coati::sampleback(work, anc, des, aln, aln.gap.len, rand);

        out << "  {\n    \"aln\": {\n";
        out << "      \"" << aln.data.names[0] << "\": ";
        out << "\"" << aln.data.seqs[0] << "\",\n";
        out << "      \"" << aln.data.names[1] << "\": ";
        out << "\"" << aln.data.seqs[1] << "\"\n";
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

}  // namespace coati
