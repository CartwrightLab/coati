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
 * @param[in] args coati::utils::args_t input parameters.
 * @param[in] aln coati::utils::alignment_t alignment information.
 */
bool marg_alignment(coati::utils::args_t& args,
                    coati::utils::alignment_t& aln) {
    coati::Matrixf P(64, 64), p_marg;
    std::ofstream out_w;

    // score alignment
    if(args.score) {
        std::cout << alignment_score(args, aln.subst_matrix) << std::endl;
        return true;
    }

    // check that length of ref sequence is multiple of 3 and gap unit size
    size_t len_a = args.fasta.seqs[0].length();
    if((len_a % 3 != 0) && (len_a % args.gap.len != 0)) {
        throw std::invalid_argument(
            "Length of reference sequence must be multiple of 3.");
    }
    // check that length of descendant sequence is multiple of gap unit size
    if(args.fasta.seqs[1].length() % args.gap.len != 0) {
        throw std::invalid_argument(
            "Length of descendant sequence must be multiple of " +
            std::to_string(args.gap.len) + ".");
    }

    // encode sequences
    auto anc = args.fasta.seqs[0];
    auto des = args.fasta.seqs[1];
    coati::utils::sequence_pair_t seq_pair =
        coati::utils::marginal_seq_encoding(anc, des);

    // dynamic programming pairwise alignment and traceback
    coati::align_pair_work_t work;
    coati::align_pair(work, seq_pair[0], seq_pair[1], aln.subst_matrix, args);
    coati::traceback(work, anc, des, aln, args.gap.len);

    if(!args.weight_file.empty()) {  // save weight and filename
        out_w.open(args.weight_file, std::ios::app | std::ios::out);
        out_w << args.fasta.path << "," << args.model << "," << aln.weight
              << std::endl;
        out_w.close();
    }

    // write alignment
    if(aln.fasta.path.extension() == ".fasta" ||
       aln.fasta.path.extension() == ".fa") {
        return coati::write_fasta(aln.fasta);
    }
    return coati::write_phylip(aln.fasta);
}

/// @private
TEST_CASE("marg_alignment") {
    SUBCASE("Alignment - output fasta") {
        coati::utils::args_t args;
        args.fasta = coati::fasta_t("", {"1", "2"}, {"CTCTGGATAGTG", "CTATAGTG"});
        args.model = "m-coati";
        args.weight_file = "score.log";
        args.output = "test-marg_alignment-fasta.fasta";

        coati::fasta_t result(args.output);
        coati::utils::alignment_t aln;
        aln.fasta.path = args.output;
        aln.fasta.names = args.fasta.names;
        coati::utils::set_subst(args, aln);

        if(std::filesystem::exists(args.output)) {
            std::filesystem::remove(args.output);
        }
        if(std::filesystem::exists(args.weight_file)) {
            std::filesystem::remove(args.weight_file);
        }

        REQUIRE(marg_alignment(args, aln));
        result = coati::read_fasta(args.output.string());

        CHECK(std::filesystem::remove(args.output));

        CHECK(result.names[0] == "1");
        CHECK(result.names[1] == "2");

        CHECK(result.seqs[0] == "CTCTGGATAGTG");
        CHECK(result.seqs[1] == "CT----ATAGTG");

        std::ifstream infile(args.weight_file);
        std::string s;
        infile >> s;
        CHECK(std::filesystem::remove(args.weight_file));
        CHECK(s.substr(s.length() - 7) == "1.51294");
    }

    SUBCASE("Alignment - output phylip") {
        coati::utils::args_t args;
        args.fasta = coati::fasta_t("", {"1", "2"}, {"GCGACTGTT", "GCGATTGCTGTT"});
        args.model = "m-coati";
        args.weight_file = "";
        args.output = "test-marg_alignment-phylip.phy";

        coati::fasta_t result(args.output);
        coati::utils::alignment_t aln;
        aln.fasta.path = args.output;
        aln.fasta.names = args.fasta.names;
        coati::utils::set_subst(args, aln);

        if(std::filesystem::exists(args.output)) {
            std::filesystem::remove(args.output);
        }

        REQUIRE(marg_alignment(args, aln));

        std::ifstream infile(args.output);
        std::string s1, s2;

        infile >> s1 >> s2;
        CHECK(s1.compare("2") == 0);
        CHECK(s2.compare("12") == 0);

        infile >> s1 >> s2;
        CHECK(s1.compare("1") == 0);
        CHECK(s2.compare("GCGA---CTGTT") == 0);

        infile >> s1 >> s2;
        CHECK(s1.compare("2") == 0);
        CHECK(s2.compare("GCGATTGCTGTT") == 0);

        CHECK(std::filesystem::remove(args.output));
    }

    SUBCASE("Alignment 2 dels - output phylip") {
        coati::utils::args_t args;
        args.fasta = coati::fasta_t("", {"1", "2"}, {"ACGTTAAGGGGT", "ACGAAT"});
        args.model = "m-coati";
        args.weight_file = "";
        args.output = "test-marg_alignment-phylip2.phy";

        coati::fasta_t result(args.output);
        coati::utils::alignment_t aln;
        aln.fasta.path = args.output;
        aln.fasta.names = args.fasta.names;
        coati::utils::set_subst(args, aln);

        if(std::filesystem::exists(args.output)) {
            std::filesystem::remove(args.output);
        }

        REQUIRE(marg_alignment(args, aln));

        std::ifstream infile(args.output);
        std::string s1, s2;

        infile >> s1 >> s2;
        CHECK(s1.compare("2") == 0);
        CHECK(s2.compare("12") == 0);

        infile >> s1 >> s2;
        CHECK(s1.compare("1") == 0);
        CHECK(s2.compare("ACGTTAAGGGGT") == 0);

        infile >> s1 >> s2;
        CHECK(s1.compare("2") == 0);
        CHECK(s2.compare("ACG--AA----T") == 0);

        CHECK(std::filesystem::remove(args.output));
    }

    SUBCASE("Alignment with gap length multiple of 3") {
        coati::utils::args_t args;
        args.fasta = coati::fasta_t("", {"1", "2"}, {"ACGTTAAGGGGT", "ACGAAT"});
        args.model = "m-coati";
        args.weight_file = "";
        args.output = "test-marg_alignment-no-frameshifts.fasta";
        args.gap.len = 3;

        coati::fasta_t result(args.output);
        coati::utils::alignment_t aln;
        aln.fasta.path = args.output;
        aln.fasta.names = args.fasta.names;
        coati::utils::set_subst(args, aln);

        if(std::filesystem::exists(args.output)) {
            std::filesystem::remove(args.output);
        }

        REQUIRE(marg_alignment(args, aln));
        result = coati::read_fasta(args.output.string());

        CHECK(std::filesystem::remove(args.output));

        CHECK(result.names[0] == "1");
        CHECK(result.names[1] == "2");

        CHECK(result.seqs[0] == "ACG---TTAAGGGGT");
        CHECK(result.seqs[1] == "ACGAAT---------");
    }

    SUBCASE("Alignment with gap length multiple of 3 - fail") {
        coati::utils::args_t args;
        args.fasta = coati::fasta_t("", {"1", "2"}, {"GCGATTGCTGT", "GCGACTGTT"});
        args.model = "m-coati";
        args.weight_file = "";
        args.output = "test-marg_alignment-no-frameshifts.fasta";
        args.gap.len = 3;

        coati::utils::alignment_t aln;
        aln.fasta.path = args.output;
        aln.fasta.names = args.fasta.names;
        coati::utils::set_subst(args, aln);
        REQUIRE_THROWS_AS(marg_alignment(args, aln), std::invalid_argument);
    }

    SUBCASE("Score alignment") {
        coati::utils::args_t args;
        args.fasta = coati::fasta_t("", {"1", "2"}, {"CTCTGGATAGTG", "CT----ATAGTG"});
        args.model = "m-coati";
        args.weight_file = "";
        args.output = "test-marg_alignment-score.fasta";
        args.score = true;

        coati::utils::alignment_t aln;
        aln.fasta.path = args.output;
        aln.fasta.names = args.fasta.names;
        coati::utils::set_subst(args, aln);
        REQUIRE(marg_alignment(args, aln));
    }

    SUBCASE("Score alignment - fail") {
        coati::utils::args_t args;
        args.fasta = coati::fasta_t("", {"1", "2"}, {"CTCTGGATAGTG", "CTATAGTG"});
        args.model = "m-coati";
        args.weight_file = "";
        args.output = "test-marg_alignment-score.fasta";
        args.score = true;
        
        coati::fasta_t result(args.output);
        coati::utils::alignment_t aln;
        aln.fasta.path = args.output;
        aln.fasta.names = args.fasta.names;
        coati::utils::set_subst(args, aln);

        REQUIRE_THROWS_AS(marg_alignment(args, aln), std::invalid_argument);
    }
}

/**
 * \brief Score alignment using marginal model.
 *
 * @param[in] args coati::utils::args_t input parameters.
 * @param[in] p_marg coati::Matrixf substitution matrix.
 *
 * \return alignment score (float).
 */
float alignment_score(coati::utils::args_t& args, coati::Matrixf& p_marg) {
    std::vector<std::string> aln = args.fasta.seqs;

    // check that both sequences have equal length
    if(aln[0].length() != aln[1].length()) {
        throw std::invalid_argument(
            "For alignment scoring both sequences must have equal length.");
    }

    // encode desc and gap-less ref sequences for subsitution matrix access
    std::string anc{aln[0]};
    boost::erase_all(anc, "-");
    coati::utils::sequence_pair_t seq_pair =
        coati::utils::marginal_seq_encoding(anc, aln[1]);

    // calculate log(1-g) log(1-e) log(g) log(e) log(pi)
    float_t no_gap = std::log1pf(-args.gap.open);
    float_t gap_stop = std::log1pf(-args.gap.extend);
    float_t gap_open = ::logf(args.gap.open);
    float_t gap_extend = ::logf(args.gap.extend);
    std::vector<coati::float_t> pi{args.pi};
    for(size_t i = 0; i < 4; i++) {
        pi[i] = ::logf(pi[i]);
    }

    float weight{0.f};
    int state{0}, ngap{0};
    for(size_t i = 0; i < aln[0].length(); i++) {
        switch(state) {
        case 0:  // subsitution
            if(aln[0][i] == '-') {
                // insertion;
                weight += gap_open;
                state = 2;
                ngap++;
            } else if(aln[1][i] == '-') {
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
            if(aln[0][i] == '-') {
                throw std::runtime_error(
                    "Insertion after deletion is not modeled.");
            } else if(aln[1][i] == '-') {
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
            if(aln[0][i] == '-') {
                // insertion_ext
                weight += gap_extend;
                ngap++;
            } else if(aln[1][i] == '-') {
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
    coati::utils::args_t args;
    coati::Matrixf P(mg94_p(0.0133, 0.2, {0.308, 0.185, 0.199, 0.308}));
    coati::Matrixf p_marg = marginal_p(P, args.pi);

    args.fasta.seqs = {"CTCTGGATAGTG", "CT----ATAGTG"};
    REQUIRE(alignment_score(args, p_marg) == doctest::Approx(1.51294f));

    args.fasta.seqs = {"CTCT--AT", "CTCTGGAT"};
    REQUIRE(alignment_score(args, p_marg) == doctest::Approx(-0.835939f));

    args.fasta.seqs = {"CTC", "CT"};
    REQUIRE_THROWS_AS(alignment_score(args, p_marg), std::invalid_argument);
}

void marg_sample(coati::utils::args_t& args,
                    coati::utils::alignment_t& aln,
                    random_t &rand) {
    coati::Matrixf P(64, 64), p_marg;

    std::ostream *pout;
    std::ofstream outfile;
    if(args.output.empty() || args.output == "-") {
        pout = &std::cout;
    } else {
        outfile.open(args.output);
        if(!outfile) {
            throw std::invalid_argument("Opening output file" +
                                        args.output.string() + "  failed.");
        }
        pout = &outfile;
    }
    std::ostream &out = *pout;

    // check that length of ref sequence is multiple of 3 and gap unit size
    size_t len_a = args.fasta.seqs[0].length();
    if((len_a % 3 != 0) && (len_a % args.gap.len != 0)) {
        throw std::invalid_argument(
            "Length of reference sequence must be multiple of 3.");
    }
    // check that length of descendant sequence is multiple of gap unit size
    if(args.fasta.seqs[1].length() % args.gap.len != 0) {
        throw std::invalid_argument(
            "Length of descendant sequence must be multiple of " +
            std::to_string(args.gap.len) + ".");
    }

    // encode sequences
    auto anc = args.fasta.seqs[0];
    auto des = args.fasta.seqs[1];
    coati::utils::sequence_pair_t seq_pair =
        coati::utils::marginal_seq_encoding(anc, des);

    // dynamic programming pairwise alignment and traceback
    coati::align_pair_work_t work;
    coati::align_pair(work, seq_pair[0], seq_pair[1], aln.subst_matrix, args);

    out << "[" << std::endl;

    for(size_t i = 0; i < args.sample_size; ++i) {
        coati::sampleback(work, anc, des, aln, args.gap.len, args.temperature, rand);

        out << "  {\n    \"aln\": {\n";
        out << "      \"" << aln.fasta.names[0] << "\": ";
        out << "\"" << aln.fasta.seqs[0] << "\",\n";
        out << "      \"" << aln.fasta.names[1] << "\": ";
        out << "\"" << aln.fasta.seqs[1] << "\"\n";
        out << "    },\n";
        out << "    \"weight\": " << ::expf(aln.weight) << ",\n";
        out << "    \"log_weight\": " << aln.weight << "\n";
        if(i < args.sample_size-1) {
            out << "  },";
        } else {
            out << "  }";
        }
        out << std::endl;
    }

    out << "]" << std::endl;
}


}  // namespace coati
