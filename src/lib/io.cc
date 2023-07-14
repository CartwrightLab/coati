/*
# Copyright (c) 2022 Juan J. Garcia Mesa <juanjosegarciamesa@gmail.com>
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

#include <coati/eigen_matrix.hpp>
#include <coati/fasta.hpp>
#include <coati/io.hpp>
#include <coati/json.hpp>
#include <coati/mg94q.tcc>
#include <coati/mutation_coati.hpp>
#include <coati/phylip.hpp>
#include <coati/utils.hpp>
#include <unsupported/Eigen/MatrixFunctions>

namespace coati::io {

/**
 * @brief Read substitution rate matrix from a CSV file
 *
 * @details Read from file a branch length and a codon substitution rate matrix.
 *  File is expected to have 3722 lines; 1 with branch length and 3721 with
 *  the following structure: codon,codon,value (e.g. AAA,AAA,0.0015).
 *
 * @param[in] file std::string path to input file.
 *
 * @retval coati::Matrixf codon substitution matrix.
 */
coati::Matrixf parse_matrix_csv(const std::string& file) {
    float br_len{NAN};
    Matrix61f Q;
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

    // fill the Q matrix (instantaneous substitution rate matrix)
    while(std::getline(input, line)) {
        std::stringstream ss(line);
        getline(ss, vec[0], ',');
        getline(ss, vec[1], ',');
        getline(ss, vec[2], ',');
        auto cod0 = coati::utils::cod64_to_61(coati::utils::cod_int(vec[0]));
        auto cod1 = coati::utils::cod64_to_61(coati::utils::cod_int(vec[1]));
        Q(cod0, cod1) = stof(vec[2]);
        count++;
    }

    input.close();

    // if file had a different number of lines that it should (61*61=3721)
    if(count != 3721) {
        throw std::invalid_argument(
            "Error reading substitution rate CSV file. Exiting!");
    }

    // normalize and exponentiate
    Matrix61f P = (Q * br_len).exp();
    P.transposeInPlace();
    return {61, 61, P.data(), P.data()+P.size()};    
}

/// @private
// GCOVR_EXCL_START
TEST_CASE("parse_matrix_csv") {
    SUBCASE("default") {
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
        for(size_t i = 0; i < 587; i++) {
            Q[mg94_indexes[i]] = mg94Q[i];
        }

        outfile << "0.0133" << std::endl;  // branch length
        std::cout.precision(10);
        for(auto i = 0; i < 61; i++) {
            for(auto j = 0; j < 61; j++) {
                outfile << codons[i] << "," << codons[j] << "," << Q[i * 61 + j]
                        << std::endl;
            }
        }

        outfile.close();
        coati::Matrix<coati::float_t> P_test(
            parse_matrix_csv("test-marg-matrix.csv"));
        for(auto i = 0; i < 61; i++) {
            for(auto j = 0; j < 61; j++) {
                CHECK_EQ(P(i, j), doctest::Approx(P_test(i, j)));
            }
        }
        REQUIRE(std::filesystem::remove("test-marg-matrix.csv"));
    }
    SUBCASE("error opening file") {
        CHECK_THROWS_AS(parse_matrix_csv(""), std::invalid_argument);
    }
    SUBCASE("too many lines") {
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
        // extra line
        outfile << codons[0] << "," << codons[0] << "," << Q[0] << std::endl;
        outfile.close();

        CHECK_THROWS_AS(parse_matrix_csv("test-marg-matrix.csv"),
                        std::invalid_argument);
        REQUIRE(std::filesystem::remove("test-marg-matrix.csv"));
    }
}
// GCOVR_EXCL_STOP

/**
 * @brief Read sequences and names in any supported format.
 *
 * @param[in] aln coati::alignment_t alignment information.
 *
 * @retval coati::data_t names and content of sequences.
 */
coati::data_t read_input(alignment_t& aln) {
    coati::data_t input_data;
    coati::file_type_t in_type;
    if(aln.data.path.empty()) {
        in_type = {"-", ".json"};
    } else {
        in_type = coati::utils::extract_file_type(aln.data.path);
    }

    // set input stream pointer
    std::istream* pin(nullptr);
    std::ifstream infile;
    if(in_type.path.empty() || in_type.path == "-") {
        pin = &std::cin;
    } else {
        infile.open(aln.data.path);
        if(!infile) {
            throw std::invalid_argument("Opening input file " +
                                        aln.data.path.string() + " failed.");
        }
        pin = &infile;
    }
    std::istream& in = *pin;

    // call reader depending on file type
    if(in_type.type_ext == ".fa" || in_type.type_ext == ".fasta") {
        input_data = read_fasta(in, aln.is_marginal());
    } else if(in_type.type_ext == ".phy") {
        input_data = read_phylip(in, aln.is_marginal());
    } else if(in_type.type_ext == ".json") {
        input_data = read_json(in, aln.is_marginal());
    } else {
        // if not supported format found
        throw std::invalid_argument("Invalid input " + aln.data.path.string() +
                                    ".");
    }
    input_data.path = aln.data.path;
    return input_data;
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
        REQUIRE(std::filesystem::remove(filename));

        CHECK_EQ(fasta.names[0], "1");
        CHECK_EQ(fasta.names[1], "2");
        CHECK_EQ(fasta.seqs[0], "CTCTGGATAGTC");
        CHECK_EQ(fasta.seqs[1], "CTATAGTC");
    }
    SUBCASE("phylip") {
        std::string filename{"test-read-input.phy"};
        outfile.open(filename);
        REQUIRE(outfile);
        outfile << "2 12" << std::endl;
        outfile << "test-sequeCTCTGGATAGTC" << std::endl;
        outfile << "2         CTCTGGATAGTC" << std::endl;
        outfile.close();

        aln.data.path = filename;
        coati::data_t phylip = read_input(aln);
        REQUIRE(std::filesystem::remove(filename));

        CHECK_EQ(phylip.names[0], "test-seque");
        CHECK_EQ(phylip.names[1], "2");
        CHECK_EQ(phylip.seqs[0], "CTCTGGATAGTC");
        CHECK_EQ(phylip.seqs[1], "CTCTGGATAGTC");
    }
    SUBCASE("json") {
        std::string filename{"test-read-input.json"};
        outfile.open(filename);
        REQUIRE(outfile);
        outfile << R"({
  "alignment": {
    "a": "CTCTGGATAGTC",
    "b": "CTATAGTC"
  },
  "score": 0.1
}
)";
        outfile.close();

        aln.data.path = filename;
        coati::data_t json = read_input(aln);
        REQUIRE(std::filesystem::remove(filename));

        CHECK_EQ(json.names[0], "a");
        CHECK_EQ(json.names[1], "b");
        CHECK_EQ(json.seqs[0], "CTCTGGATAGTC");
        CHECK_EQ(json.seqs[1], "CTATAGTC");
        CHECK_EQ(json.score, 0.1f);
    }
    SUBCASE("ext") {
        std::string filename("test-read-input.ext");
        outfile.open(filename);
        REQUIRE(outfile);
        outfile << "{\"data\":{\"names\":[\"a\",\"b\"],\"seqs\":"
                   "[\"CTCTGGATAGTC\",\"CTCTGGATAGTC\"]}}"
                << std::endl;
        outfile.close();
        aln.data.path = filename;
        CHECK_THROWS_AS(read_input(aln), std::invalid_argument);
        REQUIRE(std::filesystem::remove(filename));
    }
    SUBCASE("input file not found - fail") {
        aln.data.path = "test-read.json";
        CHECK_THROWS_AS(read_input(aln), std::invalid_argument);
    }
}
// GCOVR_EXCL_STOP

/**
 * @brief Write sequences and names in any suppported format.
 *
 * @param[in] data coati::alignment_t sequences, names, fsts, and score
 * information.
 * @param[in] aln_path coati::VectorFstStdArc FST object with alignment.
 */
void write_output(coati::alignment_t& aln) {
    coati::file_type_t out_type;
    if(aln.output.empty()) {
        out_type = {"-", ".json"};
    } else {
        out_type = coati::utils::extract_file_type(aln.output);
    }

    // set output stream pointer
    std::ostream* pout(nullptr);
    std::ofstream outfile;
    if(out_type.path == "-") {
        pout = &std::cout;
    } else {
        outfile.open(aln.output);
        pout = &outfile;
    }
    std::ostream& out = *pout;

    if(out_type.type_ext == ".fa" || out_type.type_ext == ".fasta") {
        write_fasta(aln.data, out);
    } else if(out_type.type_ext == ".phy") {
        write_phylip(aln.data, out);
    } else if(out_type.type_ext == ".json") {
        write_json(aln.data, out);
    } else {
        // not supported output format
        throw std::invalid_argument("Invalid output format " +
                                    out_type.type_ext + ".");
    }
}

/// @private
// GCOVR_EXCL_START
TEST_CASE("write_output") {
    std::vector<std::string> names = {"anc", "des"};
    std::vector<std::string> sequences = {
        "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC"
        "GTACGTACGTACGTACGTACGTACGTACGTTTTT",
        "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC"
        "GTACGTACGTACGTACGTACGTACGTACGTTTTT"};  // length > 100 to test new line

    SUBCASE("fasta") {
        coati::alignment_t aln;
        aln.data = coati::data_t("", names, sequences);
        aln.output = "test-write-output-fasta.fasta";
        write_output(aln);

        std::ifstream infile(aln.output);
        std::string s;
        infile >> s;
        CHECK_EQ(s, ">anc");
        infile >> s;
        CHECK_EQ(
            s, "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT");
        infile >> s;
        CHECK_EQ(s, "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTTTTT");
        infile >> s;
        CHECK_EQ(s, ">des");
        infile >> s;
        CHECK_EQ(
            s, "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT");
        infile >> s;
        CHECK_EQ(s, "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTTTTT");
        REQUIRE(std::filesystem::remove("test-write-output-fasta.fasta"));
    }

    SUBCASE("phylip") {
        coati::alignment_t aln;
        aln.data = coati::data_t("", names, sequences);
        aln.output = "test-write-output-phylip.phy";
        write_output(aln);

        std::ifstream infile(aln.output);
        std::string s;
        getline(infile, s);
        CHECK_EQ(s, "2 104");
        getline(infile, s);
        CHECK_EQ(
            s, "anc       ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC");
        getline(infile, s);
        CHECK_EQ(
            s, "des       ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC");
        getline(infile, s);  // empty line
        getline(infile, s);
        CHECK_EQ(s, "GTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTTTTT");
        CHECK_EQ(s, "GTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTTTTT");
        REQUIRE(std::filesystem::remove("test-write-output-phylip.phy"));
    }

    SUBCASE("json") {
        coati::alignment_t aln;
        aln.data = coati::data_t("", names, sequences);
        aln.output = "test-write-output-phylip.json";
        write_output(aln);

        std::ifstream infile(aln.output);
        std::stringstream ss;
        ss << infile.rdbuf();
        std::string s1 = ss.str();
        CHECK_EQ(s1, R"({
  "alignment": {
    "anc": "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTTTTT",
    "des": "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTTTTT"
  },
  "score": 0.0
}
)");
        REQUIRE(std::filesystem::remove("test-write-output-phylip.json"));
    }

    SUBCASE("ext") {
        coati::alignment_t aln;
        aln.data = coati::data_t("", names, sequences);
        aln.output = "test-write-output-ext.ext";
        REQUIRE_THROWS_AS(write_output(aln), std::invalid_argument);
        REQUIRE(std::filesystem::remove("test-write-output-ext.ext"));
    }
}
// GCOVR_EXCL_STOP
}  // namespace coati::io
