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
int cod_int(const std::string& codon) {
    unsigned char pos0 = codon[0];
    unsigned char pos1 = codon[1];
    unsigned char pos2 = codon[2];

    return (nt16_table[pos0] << 4) | (nt16_table[pos1] << 2) | nt16_table[pos2];
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
// GCOVR_EXCL_START
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
// GCOVR_EXCL_STOP

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
        app.add_option("reference", args.aln.ref, "Reference sequence")
            ->required();
    }
    app.add_option("-m,--model", args.aln.model, "Substitution model (coati ecm dna m-coati m-ecm)");
    if(command == Command::ALIGNPAIR || command == Command::SAMPLE) {
        app.add_option("-t,--time", args.aln.br_len,
                       "Evolutionary time/branch length")
            ->check(CLI::PositiveNumber);
    }
    if(command == Command::ALIGNPAIR) {
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
        app.add_option("-n,--sample-size", args.sample_size, "Sample size");
        app.add_option("-d, --seed", args.seeds,
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
sequence_pair_t marginal_seq_encoding(const std::string& anc,
                                      const std::string& des) {
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
    for(auto it = anc.cbegin(); it != anc.cend(); it++) {
        auto c0 = static_cast<unsigned char>(
                      nt16_table[static_cast<unsigned char>(*it)])
                  << 4;
        it++;
        auto c1 = static_cast<unsigned char>(
                      nt16_table[static_cast<unsigned char>(*it)])
                  << 2;
        it++;
        auto c2 = static_cast<unsigned char>(
            nt16_table[static_cast<unsigned char>(*it)]);
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

    CHECK(result[0][0] == static_cast<unsigned char>(0));
    CHECK(result[0][1] == static_cast<unsigned char>(1));
    CHECK(result[0][2] == static_cast<unsigned char>(2));
    CHECK(result[0][3] == static_cast<unsigned char>(126));
    CHECK(result[0][4] == static_cast<unsigned char>(127));
    CHECK(result[0][5] == static_cast<unsigned char>(128));
    CHECK(result[0][6] == static_cast<unsigned char>(189));
    CHECK(result[0][7] == static_cast<unsigned char>(190));
    CHECK(result[0][8] == static_cast<unsigned char>(191));
    CHECK(result[0][9] == static_cast<unsigned char>(63));
    CHECK(result[0][10] == static_cast<unsigned char>(64));
    CHECK(result[0][11] == static_cast<unsigned char>(65));
    CHECK(result[1][0] == static_cast<unsigned char>(0));
    CHECK(result[1][1] == static_cast<unsigned char>(1));
    CHECK(result[1][2] == static_cast<unsigned char>(2));
    CHECK(result[1][3] == static_cast<unsigned char>(3));
    CHECK(result[1][4] == static_cast<unsigned char>(4));
    CHECK(result[1][5] == static_cast<unsigned char>(5));
    CHECK(result[1][6] == static_cast<unsigned char>(6));
    CHECK(result[1][7] == static_cast<unsigned char>(7));
    CHECK(result[1][8] == static_cast<unsigned char>(8));
    CHECK(result[1][9] == static_cast<unsigned char>(9));
    CHECK(result[1][10] == static_cast<unsigned char>(10));
    CHECK(result[1][11] == static_cast<unsigned char>(11));
    CHECK(result[1][12] == static_cast<unsigned char>(12));
    CHECK(result[1][13] == static_cast<unsigned char>(13));
    CHECK(result[1][14] == static_cast<unsigned char>(14));
    CHECK(result[1][15] == static_cast<unsigned char>(15));

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
        P = parse_matrix_csv(aln.rate);
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
    auto test = [](std::string filename, file_type_t expected) {
        CAPTURE(std::move(filename));
        auto test = extract_file_type(std::move(filename));
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
 * \brief Read sequences and names for any supported format.
 *
 * @param[in] aln coati::alignment_t alignment information.
 *
 * \return coati::data_t object.
 */
data_t read_input(alignment_t& aln) {
    if(aln.output.empty()) {  // default output: json format & stdout
        aln.output = "json:-";
    }
    data_t input_data;
    coati::file_type_t in_type = coati::utils::extract_file_type(aln.data.path);

    // call reader depending on file type
    if(in_type.type_ext == ".fa" || in_type.type_ext == ".fasta") {
        input_data = read_fasta(aln.data.path, aln.is_marginal());
        input_data.out_file = coati::utils::extract_file_type(aln.output);
        return input_data;
    }
    if(in_type.type_ext == ".phy") {
        input_data = read_phylip(aln.data.path, aln.is_marginal());
        input_data.out_file = coati::utils::extract_file_type(aln.output);
        return input_data;
    }
    if(in_type.type_ext == ".json") {
        input_data = read_json(aln.data.path, aln.is_marginal());
        input_data.out_file = coati::utils::extract_file_type(aln.output);
        return input_data;
    }
    if(aln.data.path.empty()) {
        aln.data.path = "(empty)";
    }
    throw std::invalid_argument("Invalid input " + aln.data.path.string() +
                                ".");
}

/// @private
// GCOVR_EXCL_START
TEST_CASE("read_input") {
    std::ofstream outfile;
    coati::alignment_t aln;

    SUBCASE("fasta") {
        std::string filename{"test-read-input.fasta"};
        outfile.open(filename);
        REQUIRE(outfile);
        outfile << "; comment line" << std::endl;
        outfile << ">1" << std::endl << "CTCTGGATAGTC" << std::endl;
        outfile << ">2" << std::endl << "CTATAGTC" << std::endl;
        outfile.close();

        aln.data.path = filename;
        coati::data_t fasta = read_input(aln);
        CHECK(std::filesystem::remove(filename));

        CHECK(fasta.names[0] == "1");
        CHECK(fasta.names[1] == "2");
        CHECK(fasta.seqs[0] == "CTCTGGATAGTC");
        CHECK(fasta.seqs[1] == "CTATAGTC");
    }
    SUBCASE("phylip") {
        std::string filename{"test-read-input.phy"};
        outfile.open(filename);
        REQUIRE(outfile);
        outfile << "2 12" << std::endl;
        outfile << "1\tCTCTGGATAGTC" << std::endl;
        outfile << "2\tCTCTGGATAGTC" << std::endl;
        outfile.close();

        aln.data.path = filename;
        coati::data_t phylip = read_input(aln);
        CHECK(std::filesystem::remove(filename));

        CHECK(phylip.names[0] == "1");
        CHECK(phylip.names[1] == "2");
        CHECK(phylip.seqs[0] == "CTCTGGATAGTC");
        CHECK(phylip.seqs[1] == "CTCTGGATAGTC");
    }
    SUBCASE("json") {
        std::string filename{"test-read-input.json"};
        outfile.open(filename);
        REQUIRE(outfile);
        outfile << "{\"data\":{\"names\":[\"a\",\"b\"],\"seqs\":"
                   "[\"CTCTGGATAGTC\",\"CTCTGGATAGTC\"]}}"
                << std::endl;
        outfile.close();

        aln.data.path = filename;
        coati::data_t json = read_input(aln);
        CHECK(std::filesystem::remove(filename));

        CHECK(json.names[0] == "a");
        CHECK(json.names[1] == "b");
        CHECK(json.seqs[0] == "CTCTGGATAGTC");
        CHECK(json.seqs[1] == "CTCTGGATAGTC");
    }
    SUBCASE("ext") {
        std::string filename("test-read-input.ext");
        aln.data.path = filename;
        REQUIRE_THROWS_AS(read_input(aln), std::invalid_argument);
    }
}
// GCOVR_EXCL_STOP

/**
 * \brief Write sequences and names in any suppported format.
 *
 * @param[in] data coati::data_t sequences, names, fsts, and weight information.
 * @param[in] aln_path coati::VectorFstStdArc FST object with alignment.
 *
 * \return coati::data_t object.
 */
bool write_output(data_t& data, const VectorFstStdArc& aln_path) {
    // call writer depending on file type
    if(data.out_file.type_ext == ".fa" || data.out_file.type_ext == ".fasta") {
        return write_fasta(data, aln_path);
    } else if(data.out_file.type_ext == ".phy") {
        return write_phylip(data, aln_path);
    } else if(data.out_file.type_ext == ".json") {
        return write_json(data, aln_path);
    } else {
        throw std::invalid_argument("Invalid output " + data.out_file.path +
                                    ".");
    }
}

/// @private
// GCOVR_EXCL_START
TEST_CASE("write_output") {
    coati::data_t data;
    std::vector<std::string> names = {"anc", "des"};
    std::vector<std::string> sequences = {
        "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC"
        "GTACGTACGTACGTACGTACGTACGTACGTTTTT",
        "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC"
        "GTACGTACGTACGTACGTACGTACGTACGTTTTT"};  // length > 100 to test new line

    SUBCASE("fasta") {
        data = coati::data_t("", names, sequences);
        data.out_file.path = "test-write-output-fasta.fasta";
        data.out_file.type_ext = ".fasta";
        REQUIRE(write_output(data));

        std::ifstream infile(data.out_file.path);
        std::string s;
        infile >> s;
        CHECK(s == ">anc");
        infile >> s;
        CHECK(s ==
              "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT");
        infile >> s;
        CHECK(s == "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTTTTT");
        infile >> s;
        CHECK(s == ">des");
        infile >> s;
        CHECK(s ==
              "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT");
        infile >> s;
        CHECK(s == "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTTTTT");
        CHECK(std::filesystem::remove("test-write-output-fasta.fasta"));
    }

    SUBCASE("phylip") {
        data = coati::data_t("", names, sequences);
        data.out_file.path = "test-write-output-phylip.phy";
        data.out_file.type_ext = ".phy";
        REQUIRE(write_output(data));

        std::ifstream infile(data.out_file.path);
        std::string s1, s2;
        infile >> s1 >> s2;
        CHECK(s1 == "2");
        CHECK(s2 == "104");
        infile >> s1 >> s2;
        CHECK(s1 == "anc");
        CHECK(s2 ==
              "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
              "ACGTACGTACGTACGTACGTACGTA");
        infile >> s1 >> s2;
        CHECK(s1 == "des");
        CHECK(s2 ==
              "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
              "ACGTACGTACGTACGTACGTACGTA");
        // infile >> s1;  // blank line
        infile >> s1 >> s2;
        CHECK(s1 == "CGTACGTACGTTTTT");
        CHECK(s2 == "CGTACGTACGTTTTT");
        CHECK(std::filesystem::remove("test-write-output-phylip.phy"));
    }

    SUBCASE("json") {
        data = coati::data_t("", names, sequences);
        data.out_file.path = "test-write-output-phylip.json";
        data.out_file.type_ext = ".json";
        REQUIRE(write_json(data));

        std::ifstream infile(data.out_file.path);
        std::string s1;
        infile >> s1;
        CHECK(s1 ==
              "{\"data\":{\"names\":[\"anc\",\"des\"],\"seqs\":"
              "[\"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTA"
              "CGTACGTACGTACGTACGTACGTACGTACGTACGTACGTTTTT\","
              "\"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC"
              "GTACGTACGTACGTACGTACGTACGTACGTACGTACGTTTTT\"]}}");
        CHECK(std::filesystem::remove("test-write-output-phylip.json"));
    }

    SUBCASE("ext") {
        data = coati::data_t("", names, sequences);
        data.out_file.path = "test-write-output-ext.ext";
        data.out_file.type_ext = ".ext";
        REQUIRE_THROWS_AS(write_output(data), std::invalid_argument);
    }
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
    for(int i = 0; i < (aln.NumStates() - 1); siter.Next(), i++) {
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
