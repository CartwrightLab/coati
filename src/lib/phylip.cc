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

#include <coati/phylip.hpp>

namespace coati {
/**
 * @brief Read phylip format file.
 *
 * @param[in] in std::istream input stream pointing to stdin or file.
 * @param[in] marginal bool true if marginal model (if false create FSAs).
 *
 * @retval coati::data_t file path, names, and sequences.
 */
coati::data_t read_phylip(std::istream& in, bool marginal) {
    coati::data_t phylip;

    int n_seqs{0}, len_seqs{0};
    std::string line;

    // get number of sequences and length
    in >> line;
    n_seqs = std::stoi(line);
    in >> line;
    len_seqs = std::stoi(line);

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
        std::string name{line.substr(0, 10)};
        name.erase(remove_if(name.begin(), name.end(), ::isspace), name.end());
        phylip.names[i] = name;
        std::string seq{line.substr(10, line.length())};
        seq.erase(remove_if(seq.begin(), seq.end(), ::isspace), seq.end());
        phylip.seqs[i] = seq;
    }

    size_t count{0};
    // read rest of sequences
    while(in.good()) {
        size_t index = count % n_seqs;
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
                throw std::runtime_error(
                    "Creating acceptor from input phylip file failed. "
                    "Exiting!");
            }
            phylip.fsts[i] = accept;  // Add FSA
        }
    }

    return phylip;
}

/// @private
// GCOVR_EXCL_START
TEST_CASE("read_phylip") {
    std::ofstream outfile;
    std::ifstream in;
    std::string filename{"test-read-phylip.phy"};

    SUBCASE("Read test-read-phylip-fst.phy") {
        outfile.open(filename);
        REQUIRE(outfile);
        outfile << " 2 100" << std::endl;
        outfile << "VeryLongNa"
                   "CTCTGGATAGCTCTGGATAGCTCTGGATAGCTCTGGATAGCTCTGGATAGCTCTGGATA"
                   "G"
                << std::endl;
        outfile << "2         "
                   "CTATACTATACTATACTATACTATACTATACTATACTATACTATACTATACTATACTAT"
                   "A"
                << std::endl
                << std::endl;
        outfile << "CTCTGGATAGCTCTGGATAGCTCTGGATAGCTCTGGATAG" << std::endl;
        outfile << "CTATACTATACTATACTATACTATACTATACTATACTATA" << std::endl;
        outfile << std::endl;
        outfile.close();

        in.open(filename);
        REQUIRE(in);
        coati::data_t phylip = read_phylip(in, false);
        REQUIRE(std::filesystem::remove(filename));

        CHECK_EQ(phylip.names[0], "VeryLongNa");
        CHECK_EQ(phylip.names[1], "2");
        CHECK_EQ(
            phylip.seqs[0],
            "CTCTGGATAGCTCTGGATAGCTCTGGATAGCTCTGGATAGCTCTGGATAGCTCTGGATAGCTCT"
            "GGATAGCTCTGGATAGCTCTGGATAGCTCTGGATAG");
        CHECK_EQ(
            phylip.seqs[1],
            "CTATACTATACTATACTATACTATACTATACTATACTATACTATACTATACTATACTATACTAT"
            "ACTATACTATACTATACTATACTATACTATACTATA");

        CHECK_EQ(phylip.fsts[0].NumStates(), 101);
        CHECK_EQ(phylip.fsts[1].NumStates(), 101);

        for(int i = 0; i < 99; i++) {
            CHECK_EQ(phylip.fsts[0].NumArcs(i), 1);
            CHECK_EQ(phylip.fsts[1].NumArcs(i), 1);
        }
    }

    SUBCASE("Read test-read-phylip.phy") {
        outfile.open(filename);
        REQUIRE(outfile);
        outfile << "2 100" << std::endl;
        outfile << "1         "
                   "CTCTGGATAGCTCTGGATAGCTCTGGATAGCTCTGGATAGCTCTGGATAGCTCTGGATA"
                   "G"
                << std::endl;
        outfile << "2         "
                   "CTATACTATACTATACTATACTATACTATACTATACTATACTATACTATACTATACTAT"
                   "A"
                << std::endl
                << std::endl;
        outfile << "CTCTGGATAGCTCTGGATAGCTCTGGATAGCTCTGGATAG" << std::endl;
        outfile << "CTATACTATACTATACTATACTATACTATACTATACTATA" << std::endl;
        outfile << std::endl;
        outfile.close();

        in.open(filename);
        REQUIRE(in);
        coati::data_t phylip = read_phylip(in, true);
        REQUIRE(std::filesystem::remove(filename));

        CHECK_EQ(
            phylip.seqs[0],
            "CTCTGGATAGCTCTGGATAGCTCTGGATAGCTCTGGATAGCTCTGGATAGCTCTGGATAGCTCT"
            "GGATAGCTCTGGATAGCTCTGGATAGCTCTGGATAG");
        CHECK_EQ(
            phylip.seqs[1],
            "CTATACTATACTATACTATACTATACTATACTATACTATACTATACTATACTATACTATACTAT"
            "ACTATACTATACTATACTATACTATACTATACTATA");

        CHECK_EQ(phylip.names[0], "1");
        CHECK_EQ(phylip.names[1], "2");
    }
}
// GCOVR_EXCL_STOP

/**
 * @brief Write alignment in PHYLIP format.
 *
 * @param[in] phylip coati::data_t path, names, and sequences.
 * @param[in] out std::ostream output stream pointing to stdout or file.
 * @param[in] aln fst::VectorFst<fst::StdArc> aligned sequences in FST form.
 */
void write_phylip(coati::data_t& phylip, std::ostream& out,
                  const VectorFstStdArc& aln) {
    if(aln.NumStates() > 1) {
        coati::utils::fst_to_seqs(phylip, aln);
    }

    // write number of sequences and length
    out << phylip.size() << " " << phylip.seqs[0].length() << std::endl;

    // write aligned sequences to file
    size_t i = 50;
    for(size_t j = 0; j < phylip.size(); j++) {
        // write first 10 chars of name or name + spaces to fill up to 10 chars
        std::string seq_name = phylip.names[j].substr(0, 10);
        seq_name.append(10 - seq_name.length(), ' ');
        out << seq_name << phylip.seqs[j].substr(0, i) << std::endl;
    }
    out << std::endl;

    // write rest of the sequences 60 characters per line
    for(; i < phylip.seqs[0].length(); i += 60) {
        for(size_t j = 0; j < phylip.size(); j++) {
            out << phylip.seqs[j].substr(i, 60) << std::endl;
        }
        out << std::endl;
    }
}

/// @private
// GCOVR_EXCL_START
TEST_CASE("write_phylip") {
    std::ofstream outfile;
    std::string filename{"test-write-phylip.phy"};
    outfile.open(filename);
    REQUIRE(outfile);
    std::ostream& out = outfile;

    SUBCASE("Short sequences") {
        coati::data_t phylip("", {"tx_1", "taxa_2"},
                             {"CTCTGGATAGTG", "CT----ATAGTG"});

        write_phylip(phylip, out);

        std::ifstream infile(filename);
        std::string s1;

        getline(infile, s1);
        CHECK_EQ(s1, "2 12");

        getline(infile, s1);
        CHECK_EQ(s1, "tx_1      CTCTGGATAGTG");

        getline(infile, s1);
        CHECK_EQ(s1, "taxa_2    CT----ATAGTG");

        REQUIRE(std::filesystem::remove(filename));
    }

    SUBCASE("Multi-line sequences") {
        coati::data_t phylip(
            "", {"1", "2"},
            {"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
             "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
             "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
             "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"});

        write_phylip(phylip, out);

        std::ifstream infile(filename);
        std::string s1;

        getline(infile, s1);
        CHECK_EQ(s1, "2 100");
        getline(infile, s1);
        CHECK_EQ(
            s1, "1         AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
        getline(infile, s1);
        CHECK_EQ(
            s1, "2         AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
        getline(infile, s1);  // empty line
        getline(infile, s1);
        CHECK_EQ(s1, "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
        getline(infile, s1);
        CHECK_EQ(s1, "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");

        REQUIRE(std::filesystem::remove(filename));
    }
    SUBCASE("fst") {
        coati::data_t phylip("", {"1", "2"});

        VectorFstStdArc fst_write;
        fst_write.AddState();
        fst_write.SetStart(0);
        add_arc(fst_write, 0, 1, 2, 2);  // C -> C
        add_arc(fst_write, 1, 2, 4, 4);  // T -> T
        add_arc(fst_write, 2, 3, 0, 2);  // - -> C
        add_arc(fst_write, 3, 4, 1, 0);  // A -> -
        fst_write.SetFinal(4, 0.0);

        write_phylip(phylip, out, fst_write);
        std::ifstream infile(filename);
        std::string s1;

        getline(infile, s1);
        CHECK_EQ(s1, "2 4");

        getline(infile, s1);
        CHECK_EQ(s1, "1         CT-A");

        getline(infile, s1);
        CHECK_EQ(s1, "2         CTC-");

        REQUIRE(std::filesystem::remove(filename));
    }
}
// GCOVR_EXCL_STOP

}  // namespace coati
