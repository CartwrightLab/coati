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
 * \brief Hamming distance between two codons
 *
 * Number of positions in which the two codons are different.
 *
 * @param[in] cod1 uint8_t encoded codon.
 * @param[in] cod2 uint8_t encoded codon.
 *
 * \return hamming distance between two codons (int).
 */
int cod_distance(uint8_t cod1, uint8_t cod2) {
    int distance = 0;

    distance += (((cod1 & 48) >> 4) == ((cod2 & 48) >> 4) ? 0 : 1);
    distance += (((cod1 & 12) >> 2) == ((cod2 & 12) >> 2) ? 0 : 1);
    distance += ((cod1 & 3) == (cod2 & 3) ? 0 : 1);

    return distance;
}

/**
 * \brief Cast codon to position in codon list (AAA->0, AAAC->1 ... TTT->63)
 *
 * @param[in] codon std::string codon.
 *
 * \return encoded codon as its position in codon list.
 */
int cod_int(const std::string_view codon) {
    unsigned char pos0 = codon[0];
    unsigned char pos1 = codon[1];
    unsigned char pos2 = codon[2];

    return (nt16_table[pos0] << 4) | (nt16_table[pos1] << 2) | nt16_table[pos2];
}

/**
 * \brief Setup command line options
 *
 * @param[in] app CLI::App command line arguments parser from CLI11.
 * @param[in,out] args coati::args_t to store the input parameters.
 * @param[in] command std::string one of coati commands (i.e. alignpair or msa).
 *
 */
void set_cli_options(CLI::App& app, coati::args_t& args,
                     const Command& command) {
    app.add_option("input", args.aln.data.path,
                   "Input file (FASTA/PHYLIP/JSON accepted)")
        ->required();
    if(command == Command::MSA) {
        app.add_option("tree", args.aln.tree, "Newick phylogenetic tree")
            ->required()
            ->check(CLI::ExistingFile);
        app.add_option("reference", args.aln.refs, "Reference sequence")
            ->required();
    }
    app.add_option("-m,--model", args.aln.model,
                   "Substitution model (coati ecm dna m-coati m-ecm)");
    if(command == Command::ALIGNPAIR || command == Command::SAMPLE) {
        app.add_option("-t,--time", args.aln.br_len,
                       "Evolutionary time/branch length")
            ->check(CLI::PositiveNumber);
    }
    if(command == Command::ALIGNPAIR) {
        app.add_option("-r,--reference", args.aln.refs,
                       "Name of reference sequence");
        app.add_option("-R,--Reference", args.aln.refn,
                       "Position of reference sequence")
            ->check(CLI::Range(0, 1));
        app.add_option("-l,--weight", args.aln.weight_file,
                       "Write alignment score to file");
        app.add_flag("-s,--score", args.aln.score, "Score alignment");
    }
    app.add_option("-o,--output", args.aln.output, "Alignment output file");
    app.add_option("-g,--gap-open", args.aln.gap.open, "Gap opening score")
        ->check(CLI::PositiveNumber);
    app.add_option("-e,--gap-extend", args.aln.gap.extend,
                   "Gap extension score")
        ->check(CLI::PositiveNumber);
    app.add_option("-w,--omega", args.aln.omega,
                   "Nonsynonymous-synonymous bias")
        ->check(CLI::PositiveNumber);
    app.add_option("-p,--pi", args.aln.pi, "Nucleotide frequencies (A C G T)")
        ->expected(4);
    app.add_option("-k,--gap-len", args.aln.gap.len, "Set gap unit size");
    app.add_option("-x,--sigma", args.aln.sigma,
                   "GTR sigma parameters (AC AG AT CG CT GT)")
        ->expected(6);
    // specify string->value mappings
    std::map<std::string, coati::AmbiguousNucs> amb_map{
        {"AVG", coati::AmbiguousNucs::AVG},
        {"BEST", coati::AmbiguousNucs::BEST}};
    // CheckedTransformer translates and checks whether the results are either
    // in one of the strings or in one of the translations already
    app.add_option("-a,--ambiguous", args.aln.amb,
                   "Ambiguous nucleotides model", "AVG")
        ->transform(CLI::CheckedTransformer(amb_map, CLI::ignore_case));
    if(command == Command::SAMPLE) {
        // app.add_option("-T,--temperature", args.temperature, "Sampling
        // temperature");
        app.add_option("-n,--sample-size", args.sample.sample_size,
                       "Sample size");
        app.add_option("-d, --seed", args.sample.seeds,
                       "Space separated list of seed(s) used for sampling");
        //->expected(Nmin,Nmax);
    }
}

/**
 * \brief Setup command line options for coati format.
 *
 * @param[in] app CLI::App command line arguments parser from CLI11.
 * @param[in,out] args coati::args_t to store the input parameters.
 * @param[in] command std::string one of coati commands (i.e. alignpair or msa).
 *
 */
void set_cli_options_format(CLI::App& app, coati::args_t& args) {
    app.add_option("input", args.aln.data.path,
                   "Input file (FASTA/PHYLIP/JSON accepted)")
        ->required();
    app.add_option("-o,--output", args.aln.output, "Alignment output file");
    auto* phase = app.add_flag("-p,--preserve-phase",
                               args.format.preserve_phase, "Preserve phase");
    app.add_option("-c,--padding", args.format.padding,
                   "Padding char to format preserve phase")
        ->needs(phase);
    app.add_option("-s,--cut-sequences", args.format.seqs,
                   "Name of sequences to extract");
    app.add_option("-x,--cut-position", args.format.pos,
                   "Position of sequences to extract (1 based)");
}

/**
 * \brief Encode two sequences as vector<unsigned char>
 *
 * Encode ancestor (ref) sequence as codon \& phase, descendant as nucleotide.
 *  ref: AAA \& position 0 -> 0, AAA \& 1 -> 1, ... , TTT \& 2 -> 191.
 *  des: A -> 0, C -> 1, G -> 2, T -> 3.
 *
 * @param[in] anc std::string sequence of ancestor (reference).
 * @param[in] des std::string sequence of descendant.
 *
 * \return two sequences (ancestor \& descendant) encoded
 *  (coati::sequence_pair_t).
 */
sequence_pair_t marginal_seq_encoding(const std::string_view anc,
                                      const std::string_view des) {
    sequence_pair_t ret(2);
    ret[0].reserve(anc.length());
    ret[1].reserve(des.length());

    // check that REF/ancestor doesnt have ambiguous nucs
    auto pos = anc.find_first_not_of("ACGTU");
    if(pos != std::string::npos) {
        throw std::invalid_argument(
            "Ambiguous nucleotides in reference sequence not supported.");
    }

    // encode phase & codon: AAA0->0, AAA1->1, AAA2->2, AAC0->3, ... , TTT3->191
    for(size_t i = 0; i < anc.size() - 2; ++i) {
        auto c0 = static_cast<unsigned char>(
                      nt16_table[static_cast<unsigned char>(anc[i])])
                  << 4;
        i++;
        auto c1 = static_cast<unsigned char>(
                      nt16_table[static_cast<unsigned char>(anc[i])])
                  << 2;
        i++;
        auto c2 = static_cast<unsigned char>(
            nt16_table[static_cast<unsigned char>(anc[i])]);
        auto cod = (c0 | c1 | c2) * 3;
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
    std::string anc = "AAAGGGTTTCCC", des = "ACGTRYMKSWBDHVN-";
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

    anc = "AAACCCGGN";
    REQUIRE_THROWS_AS(marginal_seq_encoding(anc, des), std::invalid_argument);
    anc = "AAACCCGGR";
    REQUIRE_THROWS_AS(marginal_seq_encoding(anc, des), std::invalid_argument);
    anc = "YAACCCGGG";
    REQUIRE_THROWS_AS(marginal_seq_encoding(anc, des), std::invalid_argument);
}
// GCOVR_EXCL_STOP

/**
 * \brief Set subtitution matrix or FST according to model
 *
 * @param[in,out] aln coati::alignment_t alignment information.
 */
void set_subst(alignment_t& aln) {
    Matrixf P(64, 64);

    if(!aln.rate.empty()) {
        aln.model = "user_marg_model";
        P = coati::io::parse_matrix_csv(aln.rate);
        aln.subst_matrix = marginal_p(P, aln.pi, aln.amb);
    } else if(aln.model.compare("m-ecm") == 0) {
        P = ecm_p(aln.br_len, aln.omega);
        aln.subst_matrix = marginal_p(P, aln.pi, aln.amb);
    } else if(aln.model.compare("m-coati") == 0) {  // m-coati
        P = mg94_p(aln.br_len, aln.omega, aln.pi);
        aln.subst_matrix = marginal_p(P, aln.pi, aln.amb);
    } else if(aln.model.compare("coati") == 0) {
        aln.subst_fst = mg94(aln.br_len, aln.omega, aln.pi);
    } else if(aln.model.compare("dna") == 0) {
        aln.subst_fst = dna(aln.br_len, aln.omega, aln.pi);
    } else if(aln.model.compare("ecm") == 0) {
        aln.subst_fst = ecm(aln.br_len, aln.omega);
        aln.pi = {0.2676350, 0.2357727, 0.2539630, 0.2426323};
    } else {
        throw std::invalid_argument("Mutation model unknown.");
    }
}

/**
 * \brief Extract file type from path.
 *
 * Extracts extension and filename from both file.ext and ext:file.foo
 * Trims whitespace as well.
 *
 * @param[in] path std::string path to input file.
 *
 * \return coati::file_type_t object containing the path and extension.
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
 * \brief Convert alignment FST to std::string sequences.
 *
 * @param[in] data coati::data_t sequences, names, fst, weight information.
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

}  // namespace coati::utils
