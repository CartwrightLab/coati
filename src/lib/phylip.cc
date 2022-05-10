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

#include <coati/phylip.hpp>

constexpr std::size_t PRINT_SIZE = 100;

namespace coati {
/**
 * \brief Read phylip format file and create FSAs.
 *
 * @param[in] f_path std::string path to input file.
 * @param[in] marginal bool true if marginal model.
 *
 * \return file path, names, and sequences (data_t).
 */
coati::data_t read_phylip(const std::string& f_path, bool marginal) {
    coati::data_t phylip(f_path);

    std::istream* pin(nullptr);
    std::ifstream infile;  // input file
    coati::file_type_t in_type = coati::utils::extract_file_type(f_path);
    if(in_type.path.empty() || in_type.path == "-") {
        pin = &std::cin;  // set to stdin
        in_type.path = "-";
    } else {
        infile.open(f_path);
        if(!infile || !infile.good()) {
            throw std::invalid_argument("Opening input file " + f_path +
                                        " failed.");
        }
        pin = &infile;  // set to file
        in_type = coati::utils::extract_file_type(f_path);
    }
    std::istream& in = *pin;

    int n_seqs{0}, len_seqs{0};
    std::string line;

    // get number of sequences and length
    getline(in, line);
    std::vector<std::string> unsplitted;
    boost::trim_if(line, boost::is_any_of("\t "));
    boost::split(unsplitted, line, boost::is_any_of("\t "),
                 boost::token_compress_on);
    n_seqs = std::stoi(unsplitted[0]);
    len_seqs = std::stoi(unsplitted[1]);

    // reserve space for sequences
    phylip.names.resize(n_seqs);
    phylip.seqs.resize(n_seqs);
    phylip.fsts.resize(n_seqs);

    // read sequences name and first nucleotides
    for(int i = 0; i < n_seqs; i++) {
        phylip.seqs[i].resize(len_seqs);
        getline(in, line);
        if(line.empty()) {
            getline(in, line);
        }
        boost::split(unsplitted, line, boost::is_any_of("\t "),
                     boost::token_compress_on);
        phylip.names[i] = unsplitted[0];
        phylip.seqs[i] = unsplitted[1];
    }

    size_t count{0}, index{0};
    std::string name;
    // read rest of sequences
    while(in.good()) {
        index = count % n_seqs;
        getline(in, line);
        if(line.empty()) {
            continue;  // omit empty lines
        }
        // Remove spaces
        line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());
        phylip.seqs[index] += line;
        count++;
    }

    // if FST alignment, create FSAs
    if(!marginal) {
        for(int i = 0; i < n_seqs; i++) {
            VectorFstStdArc accept;  // create FSA with sequence
            if(!acceptor(phylip.seqs[i], accept)) {
                throw std::runtime_error("Creating acceptor from " + f_path +
                                         " failed. Exiting!");
            }
            phylip.fsts[i] = accept;  // Add FSA
        }
    }

    return phylip;
}

/// @private
// GCOVR_EXCL_START
TEST_CASE("read_phylip") {
    // cppcheck-suppress unusedVariable
    std::vector<VectorFstStdArc> fsts;
    std::ofstream outfile;

    SUBCASE("Read test-read-phylip-fst.phy") {
        std::string filename{"test-read-phylip-fst.phy"};
        outfile.open(filename);
        REQUIRE(outfile);
        outfile << " 2\t100" << std::endl << std::endl;
        outfile << "1\t"
                << "CTCTGGATAGCTCTGGATAGCTCTGGATAGCTCTGGATAGCTCTGGATAGCTCTGGATA"
                   "GCTCTGGATAGCTCTGGATAGCTCTGGATAGCTCTGG"
                << std::endl;
        outfile << "2 \t"
                << "CTATACTATACTATACTATACTATACTATACTATACTATACTATACTATACTATACTAT"
                   "ACTATACTATACTATACTATACTATACTATACTATAC"
                << std::endl;
        outfile << std::endl;  // blank line
        outfile << "ATAG" << std::endl;
        outfile << "TATA" << std::endl;
        outfile.close();

        coati::data_t phylip = read_phylip(filename, false);
        CHECK(std::filesystem::remove(phylip.path.string()));

        CHECK(phylip.seqs[0] ==
              "CTCTGGATAGCTCTGGATAGCTCTGGATAGCTCTGGATAGCTCTGGATAGCTCTGGATAGCTCT"
              "GGATAGCTCTGGATAGCTCTGGATAGCTCTGGATAG");
        CHECK(phylip.seqs[1] ==
              "CTATACTATACTATACTATACTATACTATACTATACTATACTATACTATACTATACTATACTAT"
              "ACTATACTATACTATACTATACTATACTATACTATA");

        CHECK(phylip.fsts[0].NumStates() == 101);
        CHECK(phylip.fsts[1].NumStates() == 101);

        for(int i = 0; i < 99; i++) {
            CHECK(phylip.fsts[0].NumArcs(i) == 1);
            CHECK(phylip.fsts[1].NumArcs(i) == 1);
        }
    }

    SUBCASE("Read test-read-phylip.phy") {
        std::string filename{"test-read-phylip.phy"};
        outfile.open(filename);
        REQUIRE(outfile);
        outfile << "2\t100" << std::endl;
        outfile << "1\t  "
                << "CTCTGGATAGCTCTGGATAGCTCTGGATAGCTCTGGATAGCTCTGGATAGCTCTGGATA"
                   "GCTCTGGATAGCTCTGGATAGCTCTGGATAGCTCTGG"
                << std::endl;
        outfile << "2\t"
                << "CTATACTATACTATACTATACTATACTATACTATACTATACTATACTATACTATACTAT"
                   "ACTATACTATACTATACTATACTATACTATACTATAC"
                << std::endl;
        outfile << std::endl;  // blank line
        outfile << "ATAG" << std::endl;
        outfile << "TATA" << std::endl;
        outfile.close();

        coati::data_t phylip = read_phylip(filename, true);
        CHECK(std::filesystem::remove(phylip.path.string()));

        CHECK(phylip.seqs[0] ==
              "CTCTGGATAGCTCTGGATAGCTCTGGATAGCTCTGGATAGCTCTGGATAGCTCTGGATAGCTCT"
              "GGATAGCTCTGGATAGCTCTGGATAGCTCTGGATAG");
        CHECK(phylip.seqs[1] ==
              "CTATACTATACTATACTATACTATACTATACTATACTATACTATACTATACTATACTATACTAT"
              "ACTATACTATACTATACTATACTATACTATACTATA");

        CHECK(phylip.names[0] == "1");
        CHECK(phylip.names[1] == "2");
    }

    SUBCASE("Error opening phylip") {
        REQUIRE_THROWS_AS(read_phylip("test-9999999999.phy", false),
                          std::invalid_argument);
        REQUIRE_THROWS_AS(read_phylip("test-9999999999.phy", true),
                          std::invalid_argument);
    }
}
// GCOVR_EXCL_STOP

/**
 * \brief Write alignment in PHYLIP format
 *
 * @param[in] phylip coati::data_t path, names, and sequences.
 * @param[in] aln fst::VectorFst<fst::StdArc> aligned sequences in FST form.
 */
bool write_phylip(coati::data_t& phylip, const VectorFstStdArc& aln) {
    if(aln.NumStates() > 1) {
        coati::utils::fst_to_seqs(phylip, aln);
    }
    // set output pointer
    std::ostream* pout(nullptr);
    std::ofstream outfile;
    // coati::file_type_t out_type;
    if(phylip.out_file.path == "-" || phylip.out_file.path.empty()) {
        pout = &std::cout;
    } else {
        outfile.open(phylip.out_file.path);
        // outfile.open(phylip.output);
        if(!outfile) {
            throw std::invalid_argument("Opening output file " +
                                        phylip.out_file.path + " failed.");
        }
        pout = &outfile;
        // out_type = coati::utils::extract_file_type(phylip.output.string());
    }
    std::ostream& out = *pout;

    // write aligned sequences to file
    out << phylip.size() << " " << phylip.seqs[0].length() << std::endl;
    size_t i = PRINT_SIZE - 11;
    for(size_t j = 0; j < phylip.size(); j++) {
        std::string seq_name = phylip.names[j].substr(0, 10);
        out << seq_name << std::string(10 - seq_name.length(), ' ') << " "
            << phylip.seqs[j].substr(0, i) << std::endl;
    }
    out << std::endl;

    for(; i < phylip.seqs[0].length(); i += PRINT_SIZE) {
        for(size_t j = 0; j < phylip.size(); j++) {
            out << phylip.seqs[j].substr(i, PRINT_SIZE) << std::endl;
        }
        out << std::endl;
    }

    return true;
}

/// @private
// GCOVR_EXCL_START
TEST_CASE("write_phylip") {
    SUBCASE("Short sequences") {
        coati::data_t phylip("", {"tx_1", "taxa_2"},
                             {"CTCTGGATAGTG", "CT----ATAGTG"});
        phylip.out_file = {{"test-write-phylip.phy"}, {"phy"}};

        REQUIRE(write_phylip(phylip));

        std::ifstream infile("test-write-phylip.phy");
        std::string s1, s2;

        infile >> s1 >> s2;
        CHECK(s1.compare("2") == 0);
        CHECK(s2.compare("12") == 0);

        infile >> s1 >> s2;
        CHECK(s1.compare("tx_1") == 0);
        CHECK(s2.compare("CTCTGGATAGTG") == 0);

        infile >> s1 >> s2;
        CHECK(s1.compare("taxa_2") == 0);
        CHECK(s2.compare("CT----ATAGTG") == 0);

        CHECK(std::filesystem::remove("test-write-phylip.phy"));
    }

    SUBCASE("Multi-line sequences") {
        coati::data_t phylip(
            "", {"1", "2"},
            {"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
             "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
             "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
             "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"});

        phylip.out_file = {{"test-write-phylip-long.phy"}, {"phy"}};
        REQUIRE(write_phylip(phylip));

        std::ifstream infile("test-write-phylip-long.phy");
        std::string s1, s2;

        infile >> s1 >> s2;
        CHECK(s1.compare("2") == 0);
        CHECK(s2.compare("100") == 0);

        infile >> s1 >> s2;
        CHECK(s1.compare("1") == 0);
        CHECK(s2.compare("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
                         "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA") == 0);

        infile >> s1 >> s2;
        CHECK(s1.compare("2") == 0);
        CHECK(s2.compare("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
                         "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA") == 0);
        infile >> s1 >> s2;
        CHECK(s1.compare("AAAAAAAAAAA") == 0);
        CHECK(s2.compare("AAAAAAAAAAA") == 0);

        CHECK(std::filesystem::remove("test-write-phylip-long.phy"));
    }
    SUBCASE("fst") {
        coati::data_t phylip("", {"1", "2"});
        phylip.out_file = {{"test-write-phylip.phy"}, {"phy"}};

        VectorFstStdArc fst_write;
        fst_write.AddState();
        fst_write.SetStart(0);
        add_arc(fst_write, 0, 1, 2, 2);  // C -> C
        add_arc(fst_write, 1, 2, 4, 4);  // T -> T
        add_arc(fst_write, 2, 3, 0, 2);  // - -> C
        add_arc(fst_write, 3, 4, 1, 0);  // A -> -
        fst_write.SetFinal(4, 0.0);

        REQUIRE(write_phylip(phylip, fst_write));
        std::ifstream infile("test-write-phylip.phy");
        std::string s1, s2;

        infile >> s1 >> s2;
        CHECK(s1.compare("2") == 0);
        CHECK(s2.compare("4") == 0);

        infile >> s1 >> s2;
        CHECK(s1.compare("1") == 0);
        CHECK(s2.compare("CT-A") == 0);

        infile >> s1 >> s2;
        CHECK(s1.compare("2") == 0);
        CHECK(s2.compare("CTC-") == 0);

        CHECK(std::filesystem::remove("test-write-phylip.phy"));
    }
}
// GCOVR_EXCL_STOP

}  // namespace coati
