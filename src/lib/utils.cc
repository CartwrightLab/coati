/*
# Copyright (c) 2020-2022 Juan J. Garcia Mesa <juanjosegarciamesa@gmail.com>
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

#include <boost/algorithm/string/trim.hpp>
#include <climits>
#include <coati/utils.hpp>
#include <filesystem>

namespace coati::utils {

/**
 * @brief Hamming distance between two codons.
 *
 * @details Number of positions in which the two codons are different.
 *
 * @param[in] cod1 uint8_t encoded codon.
 * @param[in] cod2 uint8_t encoded codon.
 *
 * @retval int hamming distance between two codons.
 */
int cod_distance(uint8_t cod1, uint8_t cod2) {
    int distance = 0;

    for(int i = 0; i < 3; ++i) {
        distance += (get_nuc(cod1, i) == get_nuc(cod2, i) ? 0 : 1);
    }

    return distance;
}

/**
 * @brief Get a codon's position in the codon list (AAA->0, AAAC->1 .. TTT->63).
 *
 * @details Each nucleotide is converted to its position in nucleotide list
 * (A = 0, C = 1, G = 2, T = 3). Then we apply the following:
 * Given cod = nuc0 nuc1 nuc2,
 * position = nuc0 * 2^4 + nuc1 * 2^2 + nuc2
 * Algorithm performs this in bit operations.
 *
 * @param[in] codon std::string codon.
 *
 * @retval int encoded codon as its position in codon list.
 */
int cod_int(const std::string_view codon) {
    unsigned char pos0 = codon[0];
    unsigned char pos1 = codon[1];
    unsigned char pos2 = codon[2];

    // check that REF/ancestor doesnt have ambiguous nucs
    auto pos = codon.find_first_not_of("ACGTUacgtu");
    if(pos != std::string::npos) {
        throw std::invalid_argument(
            "Ambiguous nucleotides in reference sequence not supported.");
    }

    return (nt16_table[pos0] << 4) | (nt16_table[pos1] << 2) | nt16_table[pos2];
}

/**
 * @brief Setup command line options for coati-alignpair.
 *
 * @param[in] app CLI::App command line arguments parser from CLI11.
 * @param[in,out] args coati::args_t to store the input parameters.
 */
void set_options_alignpair(CLI::App& app, coati::args_t& args) {
    app.get_formatter()->column_width(35);
    app.add_option("input", args.aln.data.path,
                   "Input file (FASTA/PHYLIP/JSON accepted)")
        ->required();
    auto* opt_m =
        app.add_option("-m,--model", args.aln.model,
                       "Substitution model (dna tri-mg tri-ecm mar-mg mar-ecm)")
            ->group("Model parameters");
    app.add_option("--sub", args.aln.rate,
                   "File with branch lengths and codon subst matrix")
        ->excludes(opt_m)
        ->group("Advanced options");
    app.add_option("-t,--time", args.aln.br_len,
                   "Evolutionary time/branch length")
        ->check(CLI::PositiveNumber)
        ->group("Model parameters");
    auto* opt_ref =
        app.add_option("-r,--ref", args.aln.refs,
                       "Name of reference sequence (default: 1st seq)")
            ->group("Advanced options");
    app.add_flag("-v,--rev-ref", args.aln.rev,
                 "Use 2nd seq as reference (default: 1st seq)")
        ->excludes(opt_ref)
        ->group("Advanced options");
    app.add_flag("-s,--score", args.aln.score,
                 "Score input alignment and exit");
    app.add_option("-o,--output", args.aln.output, "Alignment output file");
    app.add_option("-g,--gap-open", args.aln.gap.open, "Gap opening score")
        ->check(CLI::PositiveNumber)
        ->group("Model parameters");
    app.add_option("-e,--gap-extend", args.aln.gap.extend,
                   "Gap extension score")
        ->check(CLI::PositiveNumber)
        ->group("Model parameters");
    app.add_option("-w,--omega", args.aln.omega,
                   "Nonsynonymous-synonymous bias")
        ->check(CLI::PositiveNumber)
        ->group("Advanced options");
    app.add_option("-p,--pi", args.aln.pi, "Nucleotide frequencies (A C G T)")
        ->expected(4)
        ->group("Advanced options");
    app.add_option("-k,--gap-len", args.aln.gap.len, "Gap unit length")
        ->group("Model parameters");
    app.add_option("-x,--sigma", args.aln.sigma,
                   "GTR sigma parameters (AC AG AT CG CT GT)")
        ->expected(6)
        ->group("Advanced options");
    // specify string->value mappings
    std::map<std::string, coati::AmbiguousNucs> amb_map{
        {"AVG", coati::AmbiguousNucs::AVG},
        {"BEST", coati::AmbiguousNucs::BEST}};
    // CheckedTransformer translates and checks whether the results are either
    // in one of the strings or in one of the translations already
    app.add_option("-a,--ambiguous", args.aln.amb,
                   "Ambiguous nucleotides model", "AVG")
        ->transform(CLI::CheckedTransformer(amb_map, CLI::ignore_case))
        ->group("");
}

/// private
// GCOVR_EXCL_START
TEST_CASE("parse_arguments_alignpair") {
    coati::args_t args;
    CLI::App alnpair;
    coati::utils::set_options_alignpair(alnpair, args);

    std::vector<const char*> argv;
    std::vector<std::string> cli_args = {
        "alignpair", "test.fasta", "-m",   "tri-mg",  "-t",   "0.2",   "-r",
        "A",         "-s",         "-o",   "out.phy", "-g",   "0.015", "-e",
        "0.009",     "-w",         "0.21", "-p",      "0.15", "0.35",  "0.35",
        "0.15",      "-k",         "3",    "-x",      "0.1",  "0.1",   "0.1",
        "0.1",       "0.1",        "0.1",  "-a",      "AVG"};
    argv.reserve(cli_args.size() + 1);
    for(auto& arg : cli_args) {
        argv.push_back(arg.c_str());
    }
    argv.push_back(nullptr);
    alnpair.parse(static_cast<int>(argv.size() - 1), argv.data());

    CHECK_EQ(args.aln.data.path, "test.fasta");
    CHECK_EQ(args.aln.model, "tri-mg");
    CHECK_EQ(args.aln.br_len, 0.2f);
    CHECK_EQ(args.aln.refs, "A");
    CHECK(args.aln.score);
    CHECK_EQ(args.aln.output, "out.phy");
    CHECK_EQ(args.aln.gap.open, 0.015f);
    CHECK_EQ(args.aln.gap.extend, 0.009f);
    CHECK_EQ(args.aln.omega, 0.21f);
    CHECK_EQ(args.aln.pi[0], 0.15f);
    CHECK_EQ(args.aln.pi[1], 0.35f);
    CHECK_EQ(args.aln.pi[2], 0.35f);
    CHECK_EQ(args.aln.pi[3], 0.15f);
    CHECK_EQ(args.aln.gap.len, 3);
    for(size_t i = 0; i < 6; ++i) {
        CHECK_EQ(args.aln.sigma[i], 0.1f);
    }
    CHECK_EQ(args.aln.amb, coati::AmbiguousNucs::AVG);
}
// GCOVR_EXCL_STOP

/**
 * @brief Setup command line options for coati-msa.
 *
 * @param[in] app CLI::App command line arguments parser from CLI11.
 * @param[in,out] args coati::args_t to store the input parameters.
 */
void set_options_msa(CLI::App& app, coati::args_t& args) {
    app.get_formatter()->column_width(35);
    app.add_option("input", args.aln.data.path,
                   "Input file (FASTA/PHYLIP/JSON accepted)")
        ->required();
    app.add_option("tree", args.aln.tree, "Newick phylogenetic tree")
        ->required()
        ->check(CLI::ExistingFile);
    app.add_option("reference", args.aln.refs, "Name of reference sequence")
        ->required();
    app.add_option("-m,--model", args.aln.model,
                   "Substitution model (mar-mg mar-ecm)")
        ->group("Model parameters");
    app.add_option("-o,--output", args.aln.output, "Alignment output file");
    app.add_option("-g,--gap-open", args.aln.gap.open, "Gap opening score")
        ->check(CLI::PositiveNumber)
        ->group("Model parameters");
    app.add_option("-e,--gap-extend", args.aln.gap.extend,
                   "Gap extension score")
        ->check(CLI::PositiveNumber)
        ->group("Model parameters");
    app.add_option("-w,--omega", args.aln.omega,
                   "Nonsynonymous-synonymous bias")
        ->check(CLI::PositiveNumber)
        ->group("Advanced options");
    app.add_option("-p,--pi", args.aln.pi, "Nucleotide frequencies (A C G T)")
        ->expected(4)
        ->group("Advanced options");
    app.add_option("-k,--gap-len", args.aln.gap.len, "Gap unit length")
        ->group("Model parameters");
    app.add_option("-x,--sigma", args.aln.sigma,
                   "GTR sigma parameters (AC AG AT CG CT GT)")
        ->expected(6)
        ->group("Advanced options");
    // specify string->value mappings
    std::map<std::string, coati::AmbiguousNucs> amb_map{
        {"AVG", coati::AmbiguousNucs::AVG},
        {"BEST", coati::AmbiguousNucs::BEST}};
    // CheckedTransformer translates and checks whether the results are either
    // in one of the strings or in one of the translations already
    app.add_option("-a,--ambiguous", args.aln.amb,
                   "Ambiguous nucleotides model", "AVG")
        ->transform(CLI::CheckedTransformer(amb_map, CLI::ignore_case))
        ->group("");
}

/// private
// GCOVR_EXCL_START
TEST_CASE("parse_arguments_msa") {
    coati::args_t args;
    CLI::App msa;
    coati::utils::set_options_msa(msa, args);

    // create tree.newick since function checks for it to be an existing file
    std::ofstream outfile;
    outfile.open("tree.newick");
    REQUIRE(outfile);
    outfile << "((A:0.1,B:0.1):0.1);" << std::endl;
    outfile.close();

    std::vector<const char*> argv;
    std::vector<std::string> cli_args = {
        "msa",         "test.fasta", "tree.newick",  "seqA",
        "--model",     "tri-ecm",    "--output",     "out.phy",
        "--gap-open",  "0.015",      "--gap-extend", "0.009",
        "--omega",     "0.21",       "--pi",         "0.15",
        "0.35",        "0.35",       "0.15",         "--gap-len",
        "3",           "--sigma",    "0.1",          "0.1",
        "0.1",         "0.1",        "0.1",          "0.1",
        "--ambiguous", "BEST"};
    argv.reserve(cli_args.size() + 1);
    for(auto& arg : cli_args) {
        argv.push_back(arg.c_str());
    }
    argv.push_back(nullptr);
    msa.parse(static_cast<int>(argv.size() - 1), argv.data());

    CHECK_EQ(args.aln.data.path, "test.fasta");
    CHECK_EQ(args.aln.tree, "tree.newick");
    CHECK_EQ(args.aln.refs, "seqA");
    CHECK_EQ(args.aln.model, "tri-ecm");
    CHECK_EQ(args.aln.output, "out.phy");
    CHECK_EQ(args.aln.gap.open, 0.015f);
    CHECK_EQ(args.aln.gap.extend, 0.009f);
    CHECK_EQ(args.aln.omega, 0.21f);
    CHECK_EQ(args.aln.pi[0], 0.15f);
    CHECK_EQ(args.aln.pi[1], 0.35f);
    CHECK_EQ(args.aln.pi[2], 0.35f);
    CHECK_EQ(args.aln.pi[3], 0.15f);
    CHECK_EQ(args.aln.gap.len, 3);
    for(size_t i = 0; i < 6; ++i) {
        CHECK_EQ(args.aln.sigma[i], 0.1f);
    }
    CHECK_EQ(args.aln.amb, coati::AmbiguousNucs::BEST);
    REQUIRE(std::filesystem::remove("tree.newick"));
}
// GCOVR_EXCL_STOP

/**
 * @brief Setup command line options for coati-sample.
 *
 * @param[in] app CLI::App command line arguments parser from CLI11.
 * @param[in,out] args coati::args_t to store the input parameters.
 */
void set_options_sample(CLI::App& app, coati::args_t& args) {
    app.get_formatter()->column_width(35);
    app.add_option("input", args.aln.data.path,
                   "Input file (FASTA/PHYLIP/JSON accepted)")
        ->required();
    app.add_option("-t,--time", args.aln.br_len,
                   "Evolutionary time/branch length")
        ->check(CLI::PositiveNumber)
        ->group("Model parameters");
    auto* opt_m = app.add_option("-m,--model", args.aln.model,
                                 "Substitution model (mar-mg mar-ecm)")
                      ->group("Model parameters");
    app.add_option("--sub", args.aln.rate,
                   "File with branch lengths and codon subst matrix")
        ->excludes(opt_m)
        ->group("Advanced options");
    app.add_option("-o,--output", args.aln.output, "Alignment output file");
    app.add_option("-g,--gap-open", args.aln.gap.open, "Gap opening score")
        ->check(CLI::PositiveNumber)
        ->group("Model parameters");
    app.add_option("-e,--gap-extend", args.aln.gap.extend,
                   "Gap extension score")
        ->check(CLI::PositiveNumber)
        ->group("Model parameters");
    app.add_option("-w,--omega", args.aln.omega,
                   "Nonsynonymous-synonymous bias")
        ->check(CLI::PositiveNumber)
        ->group("Advanced options");
    app.add_option("-p,--pi", args.aln.pi, "Nucleotide frequencies (A C G T)")
        ->expected(4)
        ->group("Advanced options");
    app.add_option("-k,--gap-len", args.aln.gap.len, "Gap unit length")
        ->group("Model parameters");
    app.add_option("-x,--sigma", args.aln.sigma,
                   "GTR sigma parameters (AC AG AT CG CT GT)")
        ->expected(6)
        ->group("Advanced options");
    // specify string->value mappings
    std::map<std::string, coati::AmbiguousNucs> amb_map{
        {"AVG", coati::AmbiguousNucs::AVG},
        {"BEST", coati::AmbiguousNucs::BEST}};
    // CheckedTransformer translates and checks whether the results are
    // either in one of the strings or in one of the translations already
    app.add_option("-a,--ambiguous", args.aln.amb,
                   "Ambiguous nucleotides model", "AVG")
        ->transform(CLI::CheckedTransformer(amb_map, CLI::ignore_case))
        ->group("");
    // app.add_option("-T,--temperature", args.temperature, "Sampling
    // temperature");
    app.add_option("-n,--sample-size", args.sample.sample_size, "Sample size");
    app.add_option("-s, --seed", args.sample.seeds,
                   "Space separated list of seed(s) used for sampling");
}

/// private
// GCOVR_EXCL_START
TEST_CASE("parse_arguments_sample") {
    coati::args_t args;
    CLI::App sample;
    coati::utils::set_options_sample(sample, args);

    std::vector<const char*> argv;
    std::vector<std::string> cli_args = {
        "sample",      "test.fasta", "-t",           "0.2001",
        "--model",     "mar-ecm",    "--output",     "out.phy",
        "--gap-open",  "0.015",      "--gap-extend", "0.009",
        "--omega",     "0.21",       "--pi",         "0.15",
        "0.35",        "0.35",       "0.15",         "--gap-len",
        "3",           "--sigma",    "0.1",          "0.1",
        "0.1",         "0.1",        "0.1",          "0.1",
        "--ambiguous", "BEST",       "-n",           "10",
        "-s",          "42"};
    argv.reserve(cli_args.size() + 1);
    for(auto& arg : cli_args) {
        argv.push_back(arg.c_str());
    }
    argv.push_back(nullptr);
    sample.parse(static_cast<int>(argv.size() - 1), argv.data());

    CHECK_EQ(args.aln.data.path, "test.fasta");
    CHECK_EQ(args.aln.br_len, 0.2001f);
    CHECK_EQ(args.aln.model, "mar-ecm");
    CHECK_EQ(args.aln.output, "out.phy");
    CHECK_EQ(args.aln.gap.open, 0.015f);
    CHECK_EQ(args.aln.gap.extend, 0.009f);
    CHECK_EQ(args.aln.omega, 0.21f);
    CHECK_EQ(args.aln.pi[0], 0.15f);
    CHECK_EQ(args.aln.pi[1], 0.35f);
    CHECK_EQ(args.aln.pi[2], 0.35f);
    CHECK_EQ(args.aln.pi[3], 0.15f);
    CHECK_EQ(args.aln.gap.len, 3);
    for(size_t i = 0; i < 6; ++i) {
        CHECK_EQ(args.aln.sigma[i], 0.1f);
    }
    CHECK_EQ(args.aln.amb, coati::AmbiguousNucs::BEST);
    CHECK_EQ(args.sample.sample_size, 10);
    CHECK_EQ(args.sample.seeds[0], "42");
}
// GCOVR_EXCL_STOP

/**
 * @brief Setup command line options for coati format.
 *
 * @param[in] app CLI::App command line arguments parser from CLI11.
 * @param[in,out] args coati::args_t to store the input parameters.
 *
 */
void set_options_format(CLI::App& app, coati::args_t& args) {
    app.get_formatter()->column_width(35);
    app.add_option("input", args.aln.data.path,
                   "Input file (FASTA/PHYLIP/JSON accepted)")
        ->required();
    app.add_option("-o,--output", args.aln.output, "Alignment output file");
    auto* phase = app.add_flag("-p,--preserve-phase",
                               args.format.preserve_phase, "Preserve phase");
    app.add_option("-c,--padding", args.format.padding,
                   "Padding char to format preserve phase")
        ->needs(phase);
    auto* cut_seq = app.add_option("-s,--cut-seqs", args.format.names,
                                   "Name of sequences to extract");
    app.add_option("-x,--cut-pos", args.format.pos,
                   "Position of sequences to extract (1 based)")
        ->excludes(cut_seq);
}

/// private
// GCOVR_EXCL_START
TEST_CASE("parse_arguments_format") {
    coati::args_t args;
    CLI::App format;
    coati::utils::set_options_format(format, args);

    std::vector<const char*> argv;
    std::vector<std::string> cli_args = {
        "sample", "test.fasta", "--output", "out.phy", "-p",
        "-c",     "$",          "-s",       "name1",   "name2"};
    argv.reserve(cli_args.size() + 1);
    for(auto& arg : cli_args) {
        argv.push_back(arg.c_str());
    }
    argv.push_back(nullptr);
    format.parse(static_cast<int>(argv.size() - 1), argv.data());

    CHECK_EQ(args.aln.data.path, "test.fasta");
    CHECK_EQ(args.aln.output, "out.phy");
    CHECK(args.format.preserve_phase);
    CHECK_EQ(args.format.padding, "$");
    CHECK_EQ(args.format.names.size(), 2);
    CHECK_EQ(args.format.names[0], "name1");
    CHECK_EQ(args.format.names[1], "name2");
}
// GCOVR_EXCL_STOP

/**
 * @brief Encode two sequences as vector<unsigned char>.
 *
 * @details Encode ancestor (ref) sequence as codon \& phase, descendant as
 * nucleotide. ref: AAA \& position 0 -> 0, AAA \& 1 -> 1, ... , TTT \& 2 ->
 * 191. des: A -> 0, C -> 1, G -> 2, T -> 3. Ending stop codons have been
 * removed from both sequences if present. Only descendant is allowed to have
 * early stop codons (considered artifacts).
 *
 * @param[in] anc std::string sequence of ancestor (reference).
 * @param[in] des std::string sequence of descendant.
 *
 * @return coati::sequence_pair_t two sequences (ancestor \& descendant)
 * encoded.
 */
sequence_pair_t marginal_seq_encoding(const std::string_view anc,
                                      const std::string_view des) {
    sequence_pair_t ret(2);
    ret[0].reserve(anc.length());
    ret[1].reserve(des.length());

    // encode phase & codon: AAA0->0, AAA1->1, AAA2->2, AAC0->3, ... ,
    // TTT3->191
    for(size_t i = 0; i < anc.size(); i += 3) {
        auto cod = cod_int(anc.substr(i, i + 2));
        // if stop codon - throw error
        if(cod == 48 || cod == 50 || cod == 56) {
            throw std::invalid_argument("Early stop codon in ancestor.");
        }
        cod *= 3;
        ret[0].push_back(cod);
        ret[0].push_back(cod + 1);
        ret[0].push_back(cod + 2);
    }

    //  using nt16_table that converts A->0, C->1, G->2, T->3
    for(auto nuc : des) {
        ret[1].push_back(nt16_table[static_cast<unsigned char>(nuc)]);
    }

    return ret;
}

/// @private
// GCOVR_EXCL_START
TEST_CASE("marginal_seq_encoding") {
    std::string anc = "AAAGGGTTTCCCACTAGA", des = "ACGTRYMKSWBDHVN-";
    auto result = marginal_seq_encoding(anc, des);

    CHECK_EQ(result[0][0], static_cast<unsigned char>(0));
    CHECK_EQ(result[0][1], static_cast<unsigned char>(1));
    CHECK_EQ(result[0][2], static_cast<unsigned char>(2));
    CHECK_EQ(result[0][3], static_cast<unsigned char>(126));
    CHECK_EQ(result[0][4], static_cast<unsigned char>(127));
    CHECK_EQ(result[0][5], static_cast<unsigned char>(128));
    CHECK_EQ(result[0][6], static_cast<unsigned char>(189));
    CHECK_EQ(result[0][7], static_cast<unsigned char>(190));
    CHECK_EQ(result[0][8], static_cast<unsigned char>(191));
    CHECK_EQ(result[0][9], static_cast<unsigned char>(63));
    CHECK_EQ(result[0][10], static_cast<unsigned char>(64));
    CHECK_EQ(result[0][11], static_cast<unsigned char>(65));
    CHECK_EQ(result[0][12], static_cast<unsigned char>(21));
    CHECK_EQ(result[0][13], static_cast<unsigned char>(22));
    CHECK_EQ(result[0][14], static_cast<unsigned char>(23));
    CHECK_EQ(result[0][15], static_cast<unsigned char>(24));
    CHECK_EQ(result[0][16], static_cast<unsigned char>(25));
    CHECK_EQ(result[0][17], static_cast<unsigned char>(26));

    CHECK_EQ(result[1][0], static_cast<unsigned char>(0));
    CHECK_EQ(result[1][1], static_cast<unsigned char>(1));
    CHECK_EQ(result[1][2], static_cast<unsigned char>(2));
    CHECK_EQ(result[1][3], static_cast<unsigned char>(3));
    CHECK_EQ(result[1][4], static_cast<unsigned char>(4));
    CHECK_EQ(result[1][5], static_cast<unsigned char>(5));
    CHECK_EQ(result[1][6], static_cast<unsigned char>(6));
    CHECK_EQ(result[1][7], static_cast<unsigned char>(7));
    CHECK_EQ(result[1][8], static_cast<unsigned char>(8));
    CHECK_EQ(result[1][9], static_cast<unsigned char>(9));
    CHECK_EQ(result[1][10], static_cast<unsigned char>(10));
    CHECK_EQ(result[1][11], static_cast<unsigned char>(11));
    CHECK_EQ(result[1][12], static_cast<unsigned char>(12));
    CHECK_EQ(result[1][13], static_cast<unsigned char>(13));
    CHECK_EQ(result[1][14], static_cast<unsigned char>(14));
    CHECK_EQ(result[1][15], static_cast<unsigned char>(15));

    // ambiguous nucleotides in ancestor
    anc = "AAACCCGGN";
    REQUIRE_THROWS_AS(marginal_seq_encoding(anc, des), std::invalid_argument);
    anc = "AAACCCGGR";
    REQUIRE_THROWS_AS(marginal_seq_encoding(anc, des), std::invalid_argument);
    anc = "YAACCCGGG";
    REQUIRE_THROWS_AS(marginal_seq_encoding(anc, des), std::invalid_argument);
    // stop codons
    anc = "AAATAA";
    REQUIRE_THROWS_AS(marginal_seq_encoding(anc, des), std::invalid_argument);
    anc = "AAATAGGCC";
    REQUIRE_THROWS_AS(marginal_seq_encoding(anc, des), std::invalid_argument);
    anc = "TGA";
    REQUIRE_THROWS_AS(marginal_seq_encoding(anc, des), std::invalid_argument);
}
// GCOVR_EXCL_STOP

/**
 * @brief Set subtitution matrix or FST according to model.
 *
 * @param[in,out] aln coati::alignment_t alignment information containing model
 * and substitution matrix/FST.
 */
void set_subst(alignment_t& aln) {
    Matrixf P(64, 64);

    if(!aln.rate.empty()) {
        aln.model = "user_marg_model";
        P = coati::io::parse_matrix_csv(aln.rate);
        aln.subst_matrix = marginal_p(P, aln.pi, aln.amb);
    } else if(aln.model.compare("mar-ecm") == 0) {
        P = ecm_p(aln.br_len, aln.omega);
        aln.subst_matrix = marginal_p(P, aln.pi, aln.amb);
    } else if(aln.model.compare("mar-mg") == 0) {  // marginal
        P = mg94_p(aln.br_len, aln.omega, aln.pi);
        aln.subst_matrix = marginal_p(P, aln.pi, aln.amb);
    } else if(aln.model.compare("tri-mg") == 0) {
        aln.subst_fst = mg94(aln.br_len, aln.omega, aln.pi);
    } else if(aln.model.compare("dna") == 0) {
        aln.subst_fst = dna(aln.br_len, aln.omega, aln.pi);
    } else if(aln.model.compare("tri-ecm") == 0) {
        aln.subst_fst = ecm(aln.br_len, aln.omega);
        aln.pi = {0.2676350, 0.2357727, 0.2539630, 0.2426323};
    } else {
        throw std::invalid_argument("Mutation model unknown.");
    }
}

/**
 * @brief Extract file type from path.
 *
 * @details Extracts extension and filename from both file.ext and ext:file.foo.
 * Trims whitespace as well.
 *
 * @param[in] path std::string path to input file.
 *
 * @retval coati::file_type_t object containing the path and extension.
 */
file_type_t extract_file_type(std::string path) {
    constexpr auto npos = std::string::npos;

    // trim whitespace
    boost::algorithm::trim(path);

    // Format ext:path
    auto colon = path.find_first_of(':');
    if(colon != npos && colon > 1) {
        auto filepath = path.substr(colon + 1);
        auto ext = "." + path.substr(0, colon);
        return {std::move(filepath), std::move(ext)};
    }
    std::filesystem::path fpath{path};
    return {std::move(path), fpath.extension()};
}

/// @private
// GCOVR_EXCL_START
TEST_CASE("extract_file_type") {
    // NOLINTNEXTLINE(misc-unused-parameters)
    auto test = [](std::string filename, const file_type_t& expected) {
        CAPTURE(filename);
        auto test = extract_file_type(std::move(filename));
        CHECK_EQ(test.path, expected.path);
        CHECK_EQ(test.type_ext, expected.type_ext);
    };

    test("foo.bar", {"foo.bar", ".bar"});
    test("my:foo.bar", {"foo.bar", ".my"});
    test(".bar", {".bar", ""});
    test(".", {".", ""});
    test("..", {"..", ""});
    test("my:.foo.bar", {".foo.bar", ".my"});
    test(".foo.bar", {".foo.bar", ".bar"});
    test("", {"", ""});
    test(std::string{}, {{}, {}});
    test("foo:-", {"-", ".foo"});
    test("foo:bar", {"bar", ".foo"});
    test("bar:", {"", ".bar"});
    test("c:foo.bar", {"c:foo.bar", ".bar"});

    test(" \f\n\r\t\vfoo.bar \f\n\r\t\v", {"foo.bar", ".bar"});
    test(" \f\n\r\t\vmy:foo.bar \f\n\r\t\v", {"foo.bar", ".my"});
    test(" \f\n\r\t\v.bar \f\n\r\t\v", {".bar", ""});
    test(" \f\n\r\t\v", {{}, {}});
}
// GCOVR_EXCL_STOP

/**
 * @brief Convert alignment FST to std::string sequences.
 *
 * @param[in] data coati::data_t sequences, names, fst, score information.
 * @param[in] aln coati::alignment_t alignment object.
 */
void fst_to_seqs(coati::data_t& data, const VectorFstStdArc& aln) {
    fst::SymbolTable symbols;
    fill_symbol_table(symbols);

    std::string seq1, seq2;
    fst::StateIterator<fst::StdFst> siter(aln);  // FST state iterator
    for(int i = 0; i < aln.NumStates() - 1; siter.Next(), i++) {
        fst::ArcIteratorData<fst::StdArc> info;
        aln.InitArcIterator(siter.Value(), &info);
        seq1.append(symbols.Find(info.arcs[0].ilabel));
        seq2.append(symbols.Find(info.arcs[0].olabel));
    }

    data.seqs.clear();
    data.seqs.resize(2);
    data.seqs[0] = seq1;
    data.seqs[1] = seq2;

    // map all epsilons (<eps>) to gaps (-)
    boost::replace_all(data.seqs[0], "<eps>", "-");
    boost::replace_all(data.seqs[1], "<eps>", "-");
}

/// @private
// GCOVR_EXCL_START
TEST_CASE("fst_to_seqs") {
    coati::data_t data;

    VectorFstStdArc fst;
    fst.AddState();
    fst.SetStart(0);
    add_arc(fst, 0, 1, 2, 2);  // C -> C
    add_arc(fst, 1, 2, 4, 4);  // T -> T
    add_arc(fst, 2, 3, 0, 2);  // - -> C
    add_arc(fst, 3, 4, 1, 0);  // A -> -
    fst.SetFinal(4, 0.0);

    fst_to_seqs(data, fst);

    CHECK_EQ(data.seqs[0], "CT-A");
    CHECK_EQ(data.seqs[1], "CTC-");
}
// GCOVR_EXCL_STOP

/**
 * @brief Get nucleotide from codon.
 *
 * @param[in] cod uint8_t codon.
 * @param[in] pos int position in codon = {0, 1, 2}.
 *
 * @retval uint8_t nucleotide (A = 0, C = 1, G = 2, T = 3).
 *
 */
uint8_t get_nuc(uint8_t cod, int pos) {
    std::vector<uint8_t> cod_mask = {48, 12, 3};
    std::vector<uint8_t> shift = {4, 2, 0};

    return ((cod & cod_mask[pos]) >> shift[pos]);
}

/// @private
// GCOVR_EXCL_START
TEST_CASE("get_nuc") {
    std::vector<uint8_t> nucs{0, 1, 2, 3};
    for(const auto& n1 : nucs) {  // NOLINT(clang-diagnostic-unused-variable)
        // NOLINTNEXTLINE(clang-diagnostic-unused-variable)
        for(const auto& n2 : nucs) {
            // NOLINTNEXTLINE(clang-diagnostic-unused-variable)
            for(const auto& n3 : nucs) {
                CHECK_EQ(get_nuc(16 * n1 + 4 * n2 + n3, 0), n1);
                CHECK_EQ(get_nuc(16 * n1 + 4 * n2 + n3, 1), n2);
                CHECK_EQ(get_nuc(16 * n1 + 4 * n2 + n3, 2), n3);
            }
        }
    }
}
// GCOVR_EXCL_STOP

/**
 * @brief Reorder pair of input sequences so that reference is at position zero.
 *
 * @param[in,out] aln coati::alignment_t alignment data.
 */
void order_ref(coati::alignment_t& aln) {
    if(aln.data.names[0] == aln.refs) {
        // already the first sequence: do nothing
    } else if(aln.data.names[1] == aln.refs || aln.rev) {  // swap sequences
        std::swap(aln.data.names[0], aln.data.names[1]);
        std::swap(aln.data.seqs[0], aln.data.seqs[1]);
        if(aln.data.fsts.size() > 0) {
            std::swap(aln.data.fsts[0], aln.data.fsts[1]);
        }
    } else {  // aln.refs was specified and doesn't match any seq names
        throw std::invalid_argument("Name of reference sequence not found.");
    }
}

/**
 * @brief Read and validate input sequences.
 *
 * @param[in,out] aln coati::alignment_t input sequences and alignment info.
 *
 */
void process_marginal(coati::alignment_t& aln) {
    // different functions for alignpair/msa/fst/.../ ?

    if(aln.data.size() != 2) {
        throw std::invalid_argument("Exactly two sequences required.");
    }

    // set reference sequence as first sequence
    if(!aln.refs.empty() || aln.rev) {
        order_ref(aln);
    }

    // check that length of ref is multiple of 3 and gap unit length
    size_t len_a = aln.seq(0).length();
    size_t len_b = aln.seq(1).length();
    if(len_a % 3 != 0 || len_a % aln.gap.len != 0) {
        throw std::invalid_argument(
            "Length of reference sequence must be multiple of 3 and gap unit "
            "length.");
    }

    // check that length of descendant is multiple of gap unit length
    if(len_b % aln.gap.len != 0) {
        throw std::invalid_argument(
            "Length of descendant sequence must be multiple of gap unit "
            "length.");
    }

    // handle ending stop codons
    trim_end_stops(aln.data);
}

/**
 * @brief Trim end stop codons.
 *
 * @param[in,out] aln coati::alignment_t input sequences and alignment info.
 */
void trim_end_stops(coati::data_t& data) {
    for(size_t i = 0; i < data.size(); ++i) {
        std::string seq{data.seqs[i]};
        size_t len = seq.length();
        std::string last_cod{seq.substr(len - 3)};
        if(last_cod == "TAA" || last_cod == "TAG" || last_cod == "TGA") {
            data.stops.push_back(last_cod);
            data.seqs[i].erase(len - 3);
            if(data.fsts.size() > 0) {
                int end = static_cast<int>(len);
                data.fsts[i].DeleteStates({end, end - 1, end - 2});
                data.fsts[i].SetFinal(end - 3);
            }
        } else {
            data.stops.emplace_back("");
        }
    }
}

/// @private
// GCOVR_EXCL_START
TEST_CASE("trim_end_stops") {
    auto test = [](const std::vector<std::string>& raw_seqs,     // NOLINT
                   const std::vector<std::string>& exp_seqs,     // NOLINT
                   const std::vector<std::string>& exp_stops) {  // NOLINT
        coati::data_t data;
        data.seqs = raw_seqs;
        for(size_t i = 0; i < raw_seqs.size(); ++i) {
            data.names.push_back("seq_name");
        }

        trim_end_stops(data);

        CHECK_EQ(data.seqs, exp_seqs);
        CHECK_EQ(data.stops, exp_stops);
    };
    test({"AAA", "CCC"}, {"AAA", "CCC"}, {"", ""});              // no stops
    test({"AAATAA", "AAATTT"}, {"AAA", "AAATTT"}, {"TAA", ""});  // stop on ref
    test({"AAATTT", "AAATAG"}, {"AAATTT", "AAA"}, {"", "TAG"});  // stop on des
    test({"AAATGA", "AAATGA"}, {"AAA", "AAA"}, {"TGA", "TGA"});  // stop on des
    test({"AAATAA", "AAATAG"}, {"AAA", "AAA"}, {"TAA", "TAG"});  // stop on des

    auto test_tri = [](const std::vector<std::string>& raw_seqs,     // NOLINT
                       const std::vector<std::string>& exp_seqs,     // NOLINT
                       const std::vector<std::string>& exp_stops) {  // NOLINT
        coati::data_t data;
        data.seqs = raw_seqs;
        for(const auto& seq : raw_seqs) {
            data.names.push_back("seq_name");
            VectorFstStdArc fsa;  // FSA input sequence
            coati::acceptor(seq, fsa);
            data.fsts.push_back(fsa);
        }

        trim_end_stops(data);

        CHECK_EQ(data.seqs, exp_seqs);
        CHECK_EQ(data.stops, exp_stops);

        for(size_t i = 0; i < data.size(); ++i) {
            VectorFstStdArc expected;
            acceptor(data.seqs[i], expected);
            CHECK(fst::Equal(expected, data.fsts[i]));
        }
    };
    test_tri({"AAA", "CCC"}, {"AAA", "CCC"}, {"", ""});  // no stops
    // stop on reference/ancestor
    test_tri({"AAATAA", "AAATTT"}, {"AAA", "AAATTT"}, {"TAA", ""});
    // stop on descendant
    test_tri({"AAATTT", "AAATAG"}, {"AAATTT", "AAA"}, {"", "TAG"});
    // stop on both sequences
    test_tri({"AAATGA", "AAATGA"}, {"AAA", "AAA"}, {"TGA", "TGA"});
    test_tri({"AAATAA", "AAATAG"}, {"AAA", "AAA"}, {"TAA", "TAG"});
}
// GCOVR_EXCL_STOP

/**
 * @brief Restore end stop codons post pairwise alignment.
 *
 * @details Four case scenarios:
 *  (1) no end stop codons on either sequence: nothing to do.
 *  (2) only present in one sequence: add codon back and insert 3 gaps to other
        sequence.
 *  (3) present in both and same stop codon: add codons back as 3 matches.
 *  (4) present in both but different stop codon: add codons as a
 *      3 len insertion + a 3 len deletion.
 *
 * @param[in,out] data coati::data_t sequences aligned, score, and trimmed end
 * stop codons.
 */
void restore_end_stops(coati::data_t& data, const coati::gap_t& gap) {
    if(data.stops.size() != 2) {
        throw std::runtime_error("Error restoring end stop codons.");
    }

    coati::float_t gap_score = -::logf(gap.open * gap.extend * gap.extend);
    gap_score -= ::logf(1.f - gap.extend);

    if(data.stops[0] == data.stops[1]) {  // cases 1 & 3
        data.seqs[0].append(data.stops[0]);
        data.seqs[1].append(data.stops[1]);
    } else if(data.stops[0].empty()) {  // case 2 - stop in descendant
        data.seqs[0].append("---");
        data.seqs[1].append(data.stops[1]);
        data.score += gap_score;
    } else if(data.stops[1].empty()) {  // case 2 - stop in ancestor
        data.seqs[0].append(data.stops[0]);
        data.seqs[1].append("---");
        data.score += gap_score;
    } else {  // case 4
        data.seqs[0].append("---" + data.stops[0]);
        data.seqs[1].append(data.stops[1] + "---");
        data.score += gap_score + gap_score;
    }
}

/// @private
// GCOVR_EXCL_START
TEST_CASE("restore_end_stops") {
    auto test = [](const std::vector<std::string>& seqs,
                   const std::vector<std::string>& stops,
                   const std::vector<std::string>& exp_seqs) {  // NOLINT
        coati::data_t data;
        coati::gap_t gap;
        data.seqs = seqs;
        data.stops = stops;

        restore_end_stops(data, gap);

        CHECK_EQ(data.seqs, exp_seqs);
    };
    test({"AAA", "AAA"}, {"TAA", "TAA"}, {"AAATAA", "AAATAA"});  // same end cod
    test({"", ""}, {"TAA", "TAA"}, {"TAA", "TAA"});              // same end cod
    test({"CGA", "CGA"}, {"", ""}, {"CGA", "CGA"});              // no stop cod
    test({"CTA", "CTA"}, {"TAG", "TGA"},
         {"CTA---TAG", "CTATGA---"});                         // diff stop cod
    test({"TGC", "TGC"}, {"", "TAA"}, {"TGC---", "TGCTAA"});  // one stop cod
    test({"TGC---", "TGCCAC"}, {"", "TAA"},
         {"TGC------", "TGCCACTAA"});                         // one stop cod
    test({"CGG", "CGG"}, {"TAG", ""}, {"CGGTAG", "CGG---"});  // one stop cod
    // data.stops size != 2
    coati::data_t d;
    coati::gap_t gap;  // NOLINT
    d.stops = {""};
    CHECK_THROWS_AS(restore_end_stops(d, gap), std::runtime_error);
}
// GCOVR_EXCL_STOP

/**
 * @brief Read and validate input sequences.
 *
 * @param[in,out] aln coati::alignment_t input sequences and alignment info.
 */
void process_triplet(coati::alignment_t& aln) {
    if(aln.data.size() != 2) {
        throw std::invalid_argument("Exactly two sequences required.");
    }

    // set reference sequence as first sequence
    if(!aln.refs.empty() || aln.rev) {
        order_ref(aln);
    }

    if(aln.seq(0).length() % 3 != 0) {
        throw std::invalid_argument(
            "Length of reference sequence must be multiple of 3.");
    }

    // no early stop codons allowed in reference/ancestor
    std::string cod;
    for(size_t i = 0; i < aln.seq(0).length() - 3; i += 3) {
        cod = aln.seq(0).substr(i, 3);
        if(cod == "TAA" || cod == "TAG" || cod == "TGA") {
            throw std::invalid_argument("Early stop codon in ancestor.");
        }
    }

    // handle ending stop codons
    trim_end_stops(aln.data);
}

}  // namespace coati::utils
