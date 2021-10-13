/*
# Copyright (c) 2020-2021 Juan J. Garcia Mesa <juanjosegarciamesa@gmail.com>
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

#include <climits>
#include <coati/utils.hpp>

#include <boost/algorithm/string/trim.hpp>
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
int cod_int(const std::string& codon) {
    unsigned char pos0 = codon[0];
    unsigned char pos1 = codon[1];
    unsigned char pos2 = codon[2];

    return (nt4_table[pos0] << 4) | (nt4_table[pos1] << 2) | nt4_table[pos2];
}

/**
 * \brief Read substitution rate matrix from a CSV file
 *
 * Read from file a branch length and a codon substitution rate matrix.
 *  File is expected to have 4067 lines; 1 with branch length and 4096 with
 *  the following structure: codon,codon,value (e.g. AAA,AAA,0.0015).
 *
 * @param[in] file std::string path to input file.
 *
 * \return codon substitution matrix (coati::Matrixf).
 */
coati::Matrixf parse_matrix_csv(const std::string& file) {
    float br_len{NAN};
    Matrix64f Q;
    std::ifstream input(file);
    if(!input.good()) {
        throw std::invalid_argument("Error opening file " + file + ".");
    }

    std::string line;
    // Read branch length
    getline(input, line);
    br_len = stof(line);

    std::vector<std::string> vec{"", "", ""};
    int count = 0;

    while(std::getline(input, line)) {
        std::stringstream ss(line);
        getline(ss, vec[0], ',');
        getline(ss, vec[1], ',');
        getline(ss, vec[2], ',');
        Q(cod_int(vec[0]), cod_int(vec[1])) = stof(vec[2]);
        count++;
    }

    input.close();

    if(count != 4096) {
        throw std::invalid_argument(
            "Error reading substitution rate CSV file. Exiting!");
    }

    Q = Q * br_len;
    Q = Q.exp();

    coati::Matrix<coati::float_t> P(64, 64, Q);

    return P;
}

/// @private
TEST_CASE("parse_matrix_csv") {
    std::ofstream outfile;
    coati::Matrix<coati::float_t> P(
        mg94_p(0.0133, 0.2, {0.308, 0.185, 0.199, 0.308}));

    const std::vector<std::string> codons = {
        "AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC",
        "AGG", "AGT", "ATA", "ATC", "ATG", "ATT", "CAA", "CAC", "CAG", "CAT",
        "CCA", "CCC", "CCG", "CCT", "CGA", "CGC", "CGG", "CGT", "CTA", "CTC",
        "CTG", "CTT", "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT",
        "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT", "TAA", "TAC",
        "TAG", "TAT", "TCA", "TCC", "TCG", "TCT", "TGA", "TGC", "TGG", "TGT",
        "TTA", "TTC", "TTG", "TTT"};

    outfile.open("test-marg-matrix.csv");
    REQUIRE(outfile);

    float Q[4096]{0.0f};
    for(auto i = 0; i < 640; i++) {
        Q[mg94_indexes[i]] = mg94Q[i];
    }

    outfile << "0.0133" << std::endl;  // branch length
    for(auto i = 0; i < 64; i++) {
        for(auto j = 0; j < 64; j++) {
            outfile << codons[i] << "," << codons[j] << "," << Q[i * 64 + j]
                    << std::endl;
        }
    }

    outfile.close();
    coati::Matrix<coati::float_t> P_test(
        parse_matrix_csv("test-marg-matrix.csv"));
    for(auto i = 0; i < 64; i++) {
        for(auto j = 0; j < 64; j++) {
            CHECK(P(i, j) == doctest::Approx(P_test(i, j)));
        }
    }
    CHECK(std::filesystem::remove("test-marg-matrix.csv"));
}

/**
 * \brief Setup command line options
 *
 * @param[in] app CLI::App command line arguments parser from CLI11.
 * @param[in,out] args coati::utils::args_t to store the input parameters.
 * @param[in] command std::string one of coati commands (i.e. alignpair or msa).
 *
 */
void set_cli_options(CLI::App& app, coati::utils::args_t& args,
                     const std::string& command) {
    // Add new options/flags
    app.add_option("fasta", args.fasta.path, "Fasta file path")
        ->required()
        ->check(CLI::ExistingFile);
    if(command == "msa") {
        app.add_option("tree", args.tree, "Newick phylogenetic tree")
            ->required()
            ->check(CLI::ExistingFile);
        app.add_option("reference", args.ref, "Reference sequence")->required();
    }
    app.add_option("-m,--model", args.model, "Substitution model");
    if(command == "alignpair" || command == "sample") {
        app.add_option("-t,--time", args.br_len,
                       "Evolutionary time/branch length")
            ->check(CLI::PositiveNumber);
    }
    if(command == "alignpair") {
        app.add_option("-l,--weight", args.weight_file,
                       "Write alignment score to file");
        app.add_flag("-s,--score", args.score, "Score alignment");
    }
    app.add_option("-o,--output", args.output, "Alignment output file");
    app.add_option("-g,--gap-open", args.gap.open, "Gap opening score")
        ->check(CLI::PositiveNumber);
    app.add_option("-e,--gap-extend", args.gap.extend, "Gap extension score")
        ->check(CLI::PositiveNumber);
    app.add_option("-w,--omega", args.omega, "Nonsynonymous-synonymous bias")
        ->check(CLI::PositiveNumber);
    app.add_option("-p,--pi", args.pi, "Nucleotide frequencies (A C G T)")
        ->expected(4);
    app.add_option("-n,--gap-len", args.gap.len, "Set gap unit size");
}

/**
 * \brief Encode two sequences as vector<unsigned char>
 *
 * Encode ancestor (ref) sequence as codon \& phase, descendant as nucleotide.
 *  ref: AAA \& position 0 -> 0, AAA \& 1 -> 1, ... , TTT \& 3 -> 191.
 *  des: A -> 0, C -> 1, G -> 2, T -> 3.
 *
 * @param[in] anc std::string sequence of ancestor (reference).
 * @param[in] des std::string sequence of descendant.
 *
 * \return two sequences (ancestor \& descendant) encoded
 *  (coati::sequence_pair_t).
 */
sequence_pair_t marginal_seq_encoding(const std::string& anc,
                                      const std::string& des) {
    sequence_pair_t ret(2);
    ret[0].reserve(anc.length());
    ret[1].reserve(des.length());

    // encode phase & codon: AAA0->0, AAA1->1, AAA2->2, AAC0->3, ... , TTT3->191
    for(auto it = anc.cbegin(); it != anc.cend(); it++) {
        auto c0 = static_cast<unsigned char>(
                      nt4_table[static_cast<unsigned char>(*it)])
                  << 4;
        it++;
        auto c1 = static_cast<unsigned char>(
                      nt4_table[static_cast<unsigned char>(*it)])
                  << 2;
        it++;
        auto c2 = static_cast<unsigned char>(
            nt4_table[static_cast<unsigned char>(*it)]);
        auto cod = (c0 | c1 | c2) * 3;
        ret[0].push_back(cod);
        ret[0].push_back(cod + 1);
        ret[0].push_back(cod + 2);
    }

    //  using nt4_table that converts A->0, C->1, G->2, T->3
    for(auto nuc : des) {
        ret[1].push_back(nt4_table[static_cast<unsigned char>(nuc)]);
    }

    return ret;
}

/// @private
TEST_CASE("marginal_seq_encoding") {
    std::string anc = "AAAGGGTTT", des = "ACGT-";
    auto result = marginal_seq_encoding(anc, des);

    CHECK(result[0][0] == static_cast<unsigned char>(0));
    CHECK(result[0][1] == static_cast<unsigned char>(1));
    CHECK(result[0][2] == static_cast<unsigned char>(2));
    CHECK(result[0][3] == static_cast<unsigned char>(126));
    CHECK(result[0][4] == static_cast<unsigned char>(127));
    CHECK(result[0][5] == static_cast<unsigned char>(128));
    CHECK(result[0][6] == static_cast<unsigned char>(189));
    CHECK(result[0][7] == static_cast<unsigned char>(190));
    CHECK(result[0][8] == static_cast<unsigned char>(191));
    CHECK(result[1][0] == static_cast<unsigned char>(0));
    CHECK(result[1][1] == static_cast<unsigned char>(1));
    CHECK(result[1][2] == static_cast<unsigned char>(2));
    CHECK(result[1][3] == static_cast<unsigned char>(3));
    CHECK(result[1][4] == static_cast<unsigned char>(4));
}

/**
 * \brief Set subtitution matrix or FST according to model
 *
 * @param[in,out] args coati::utils::args_t input arguments.
 * @param[in,out] aln coati::alignment_t alignment information.
 */
void set_subst(args_t& args, alignment_t& aln) {
    Matrixf P(64, 64);

    if(!args.rate.empty()) {
        aln.model = "user_marg_model";
        P = parse_matrix_csv(args.rate);
        aln.subst_matrix = marginal_p(P, args.pi);
    } else if(args.model.compare("m-ecm") == 0) {
        P = ecm_p(args.br_len, args.omega);
        aln.subst_matrix = marginal_p(P, args.pi);
    } else if(args.model.compare("m-coati") == 0) {  // m-coati
        P = mg94_p(args.br_len, args.omega, args.pi);
        aln.subst_matrix = marginal_p(P, args.pi);
    } else if(args.model.compare("coati") == 0) {
        aln.subst_fst = mg94(args.br_len, args.omega, args.pi);
    } else if(args.model.compare("dna") == 0) {
        aln.subst_fst = dna(args.br_len, args.omega, args.pi);
    } else if(args.model.compare("ecm") == 0) {
        aln.subst_fst = ecm(args.br_len, args.omega);
        args.pi = {0.2676350, 0.2357727, 0.2539630, 0.2426323};
    } else {
        throw std::invalid_argument("Mutation model unknown.");
    }

    aln.model = args.model;
}

// extracts extension and filename from both file.ext and ext:file.foo
// trims whitespace as well
file_type_t extract_file_type(std::string path) {
    constexpr auto npos = std::string::npos;

    // trim whitespace
    boost::algorithm::trim(path);

    // Format ext:path
    auto colon = path.find_first_of(':');
    if(colon != npos && colon > 1) {
        auto filepath = path.substr(colon+1);
        auto ext = "." + path.substr(0,colon);
        return {std::move(filepath), std::move(ext)};
    }
    std::filesystem::path fpath{path};
    return {std::move(path), fpath.extension()};
}

/// @private
TEST_CASE("extract_file_type") {
    auto test = [](std::string filename, file_type_t expected) {
        CAPTURE(filename);
        auto test = extract_file_type(filename);
        CHECK(test.path == expected.path);
        CHECK(test.type_ext == expected.type_ext);
    };

    test("foo.bar", {"foo.bar", ".bar"});
    test("my:foo.bar", {"foo.bar", ".my"});
    test(".bar", {".bar", ""});
    test(".", {".", ""});
    test("..", {"..", ""});
    test("my:.foo.bar", {".foo.bar", ".my"});
    test(".foo.bar", {".foo.bar", ".bar"});
    test("", {"", ""});
    test(std::string{}, {{},{}});
    test("foo:-", {"-", ".foo"});
    test("foo:bar", {"bar", ".foo"});
    test("bar:", {"", ".bar"});
    test("c:foo.bar", {"c:foo.bar", ".bar"});

    test(" \f\n\r\t\vfoo.bar \f\n\r\t\v", {"foo.bar", ".bar"});
    test(" \f\n\r\t\vmy:foo.bar \f\n\r\t\v", {"foo.bar", ".my"});
    test(" \f\n\r\t\v.bar \f\n\r\t\v", {".bar", ""});
    test(" \f\n\r\t\v", {{},{}});
}

}  // namespace coati::utils
