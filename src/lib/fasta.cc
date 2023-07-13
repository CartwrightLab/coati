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

#include <coati/fasta.hpp>
#include <coati/mutation_fst.hpp>

namespace coati {

/**
 * @brief Read fasta format file.
 *
 * @param[in] in std::istream input stream pointing to stdin or file.
 * @param[in] marginal bool true if model is marginal.
 *
 * @retval coati::data_t object with names and content of sequences.
 */
coati::data_t read_fasta(std::istream& in, bool marginal) {
    coati::data_t fasta;

    std::string line, name, content;
    while(in.good()) {
        getline(in, line);
        if(line.empty()) {
            continue;  // omit empty lines
        }
        if(line[0] == ';') {  // omit comment lines
            continue;
        }
        if(line[0] == '>') {  // Identifier marker
            if(!name.empty()) {
                fasta.seqs.push_back(content);
                name.clear();
            }
            // Add name of sequence
            name = line.substr(1);
            if(name.empty()) {
                throw std::invalid_argument(
                    "Input fasta file contains a sequence without a name.");
            }
            fasta.names.push_back(name);
            content.clear();
        } else if(!name.empty()) {
            // Remove spaces
            line.erase(remove_if(line.begin(), line.end(), ::isspace),
                       line.end());
            content += line;
        }
    }
    if(!name.empty()) {  // Add last sequence if needed
        fasta.seqs.push_back(content);
    }

    if(!marginal) {  // if model is not marginal, create FSA.
        for(size_t i = 0; i < fasta.seqs.size(); i++) {
            VectorFstStdArc accept;  // create FSA with sequence
            if(!acceptor(fasta.seqs[i], accept)) {
                throw std::runtime_error(
                    "Creating acceptor from input fasta file failed. Exiting!");
            }
            fasta.fsts.push_back(accept);  // Add FSA
        }
    }

    return fasta;
}

/// @private
// GCOVR_EXCL_START
TEST_CASE("read_fasta") {
    std::ofstream outfile;
    std::ifstream in;
    std::string filename{"test-read-fasta.fasta"};

    SUBCASE("Read test-read-fasta.fasta") {
        outfile.open(filename);
        REQUIRE(outfile);
        outfile << "; comment line" << std::endl;
        outfile << ">1" << std::endl << "CTCTGGATAGTC" << std::endl;
        outfile << ">2" << std::endl << "CTATAGTC" << std::endl;
        outfile.close();

        in.open(filename);
        REQUIRE(in);
        coati::data_t fasta = read_fasta(in, true);
        REQUIRE(std::filesystem::remove(filename));

        CHECK_EQ(fasta.names[0], "1");
        CHECK_EQ(fasta.names[1], "2");
        CHECK_EQ(fasta.seqs[0], "CTCTGGATAGTC");
        CHECK_EQ(fasta.seqs[1], "CTATAGTC");
    }

    SUBCASE("Read test-read-fasta-fst.fasta") {
        outfile.open(filename);
        REQUIRE(outfile);
        outfile << "; comment line" << std::endl;
        outfile << ">1" << std::endl << "CTCTGGATAGTG" << std::endl;
        outfile << ">2" << std::endl << "CTATAGTG" << std::endl;
        outfile.close();

        in.open(filename);
        REQUIRE(in);
        coati::data_t fasta = read_fasta(in, false);
        REQUIRE(std::filesystem::remove(filename));

        CHECK_EQ(fasta.seqs[0], "CTCTGGATAGTG");
        CHECK_EQ(fasta.seqs[1], "CTATAGTG");

        CHECK_EQ(fasta.fsts[0].NumStates(), 13);
        CHECK_EQ(fasta.fsts[1].NumStates(), 9);

        for(int i = 0; i < 12; i++) {
            CHECK_EQ(fasta.fsts[0].NumArcs(i), 1);
        }

        for(int i = 0; i < 8; i++) {
            CHECK_EQ(fasta.fsts[0].NumArcs(i), 1);
        }
    }

    SUBCASE("sequence before first name") {
        outfile.open(filename);
        REQUIRE(outfile);
        outfile << "CTCTGGATA" << std::endl;
        outfile << ">1" << std::endl << "CTCTGGATAGTC" << std::endl;
        outfile << ">2" << std::endl << "CTATAGTC" << std::endl;
        outfile.close();

        in.open(filename);
        REQUIRE(in);
        coati::data_t fasta = read_fasta(in, true);
        REQUIRE(std::filesystem::remove(filename));

        CHECK_EQ(fasta.names[0], "1");
        CHECK_EQ(fasta.names[1], "2");
        CHECK_EQ(fasta.seqs[0], "CTCTGGATAGTC");
        CHECK_EQ(fasta.seqs[1], "CTATAGTC");
    }
    SUBCASE("empty name - fail") {
        outfile.open(filename);
        REQUIRE(outfile);
        outfile << ">" << std::endl << "CTCTGGATAGTC" << std::endl;
        outfile << ">2" << std::endl << "CTATAGTC" << std::endl;
        outfile.close();

        in.open(filename);
        REQUIRE(in);
        CHECK_THROWS_AS(read_fasta(in, true), std::invalid_argument);
        REQUIRE(std::filesystem::remove(filename));
    }
}
// GCOVR_EXCL_STOP

/**
 * @brief Write alignment in fasta format.
 *
 * @param[in] fasta coati::data_t output path, names, and aligned sequences.
 * @param[in] out std::ostream output stream pointing to stdout or file.
 * @param[in] aln VectorFstStdArc FST with alignment (optional).
 */
void write_fasta(coati::data_t& fasta, std::ostream& out) {
    // for each aligned sequence write name and content (60 characters/line)
    for(size_t i = 0; i < fasta.size(); i++) {
        out << ">" << fasta.names[i] << std::endl;
        for(size_t j = 0; j < fasta.seqs[i].size(); j += 60) {
            out << fasta.seqs[i].substr(j, 60) << std::endl;
        }
    }
}

/// @private
// GCOVR_EXCL_START
TEST_CASE("write_fasta") {
    std::ofstream outfile;
    std::string filename{"test-write-fasta.fasta"};
    outfile.open(filename);
    REQUIRE(outfile);
    std::ostream& out = outfile;

    SUBCASE("marginal") {
        coati::data_t fasta("", {"1", "2"}, {"CTCTGGATAGTG", "CTATAGTG"});

        write_fasta(fasta, out);

        std::ifstream infile(filename);
        std::string s1;
        infile >> s1;
        CHECK_EQ(s1.compare(">1"), 0);
        infile >> s1;
        CHECK_EQ(s1.compare("CTCTGGATAGTG"), 0);
        infile >> s1;
        CHECK_EQ(s1.compare(">2"), 0);
        infile >> s1;
        CHECK_EQ(s1.compare("CTATAGTG"), 0);
        REQUIRE(std::filesystem::remove(filename));
    }
    SUBCASE("data_t size fails - diff number of names and seqs") {
        coati::data_t fasta("", {"1", "2", "3"}, {"CTC", "CTA"});
        CHECK_THROWS_AS(write_fasta(fasta, out), std::invalid_argument);
        REQUIRE(std::filesystem::remove(filename));
    }
}
// GCOVR_EXCL_STOP
}  // namespace coati
