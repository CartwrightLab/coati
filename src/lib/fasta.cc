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

#include <coati/fasta.hpp>

namespace coati {

/**
 * \brief Read fasta format file.
 *
 * @param[in] f_path std::string path to fasta file.
 * @param[in] marginal bool true if model is marginal.
 *
 * \return coati::data_t object with fasta information.
 */
coati::data_t read_fasta(const std::string& f_path, bool marginal) {
    coati::data_t fasta(f_path);

    // set input pointer and file type
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
            fasta.names.push_back(name);
            content.clear();
        } else if(!name.empty()) {
            // Remove spaces
            line.erase(remove_if(line.begin(), line.end(), ::isspace),
                       line.end());
            content += line;
        }
    }
    if(!name.empty()) {  // Add last sequence FSA if needed
        fasta.seqs.push_back(content);
    }

    if(!marginal) {  // if model is not marginal, create FSA.
        for(size_t i = 0; i < fasta.seqs.size(); i++) {
            VectorFstStdArc accept;  // create FSA with sequence
            if(!acceptor(fasta.seqs[i], accept)) {
                throw std::runtime_error("Creating acceptor from " + f_path +
                                         " failed. Exiting!");
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

    SUBCASE("Read test-read-fasta.fasta") {
        std::string filename{"test-read-fasta.fasta"};
        outfile.open(filename);
        REQUIRE(outfile);
        outfile << "; comment line" << std::endl;
        outfile << ">1" << std::endl << "CTCTGGATAGTC" << std::endl;
        outfile << ">2" << std::endl << "CTATAGTC" << std::endl;
        outfile.close();

        coati::data_t fasta = read_fasta(filename, true);
        CHECK(std::filesystem::remove(fasta.path));

        CHECK(fasta.names[0] == "1");
        CHECK(fasta.names[1] == "2");
        CHECK(fasta.seqs[0] == "CTCTGGATAGTC");
        CHECK(fasta.seqs[1] == "CTATAGTC");
    }

    SUBCASE("Read test-read-fasta-fst.fasta") {
        std::string filename{"test-read-fasta-fst.fasta"};
        outfile.open(filename);
        REQUIRE(outfile);
        outfile << "; comment line" << std::endl;
        outfile << ">1" << std::endl << "CTCTGGATAGTG" << std::endl;
        outfile << ">2" << std::endl << "CTATAGTG" << std::endl;
        outfile.close();

        coati::data_t fasta = read_fasta(filename, false);
        CHECK(std::filesystem::remove(fasta.path.string()));

        CHECK(fasta.seqs[0] == "CTCTGGATAGTG");
        CHECK(fasta.seqs[1] == "CTATAGTG");

        CHECK(fasta.fsts[0].NumStates() == 13);
        CHECK(fasta.fsts[1].NumStates() == 9);

        for(int i = 0; i < 12; i++) {
            CHECK(fasta.fsts[0].NumArcs(i) == 1);
        }

        for(int i = 0; i < 8; i++) {
            CHECK(fasta.fsts[0].NumArcs(i) == 1);
        }
    }

    SUBCASE("Error opening fasta") {
        REQUIRE_THROWS_AS(read_fasta("test-9999999999.fasta", false),
                          std::invalid_argument);
        REQUIRE_THROWS_AS(read_fasta("test-9999999999.fasta", true),
                          std::invalid_argument);
    }
}
// GCOVR_EXCL_STOP

/**
 * \brief Write alignment in fasta format.
 *
 * @param[in] fasta data_t output path, names, and aligned sequences.
 * @param[in] aln VectorFstStdArc FST with alignment (optional).
 */
bool write_fasta(coati::data_t& fasta, const VectorFstStdArc& aln) {
    if(aln.NumStates() > 0) {  // if FST alignment
        coati::utils::fst_to_seqs(fasta, aln);
    }
    // set output pointer
    std::ostream* pout(nullptr);
    std::ofstream outfile;
    if(fasta.out_file.path == "-" || fasta.out_file.path.empty()) {
        pout = &std::cout;
    } else {
        outfile.open(fasta.out_file.path);
        if(!outfile) {
            throw std::invalid_argument("Opening output file " +
                                        fasta.out_file.path + " failed.");
        }
        pout = &outfile;
    }
    std::ostream& out = *pout;

    for(size_t i = 0; i < fasta.size(); i++) {
        out << ">" << fasta.names[i] << std::endl;
        for(size_t j = 0; j < fasta.seqs[i].size(); j += 60) {
            out << fasta.seqs[i].substr(j, 60) << std::endl;
        }
    }

    return true;
}

/// @private
// GCOVR_EXCL_START
TEST_CASE("write_fasta") {
    SUBCASE("marginal") {
        coati::data_t fasta("", {"1", "2"}, {"CTCTGGATAGTG", "CTATAGTG"});
        fasta.out_file.path = "test-write-fasta.fasta";

        REQUIRE(write_fasta(fasta));

        std::ifstream infile("test-write-fasta.fasta");
        std::string s1;
        infile >> s1;
        CHECK(s1.compare(">1") == 0);
        infile >> s1;
        CHECK(s1.compare("CTCTGGATAGTG") == 0);
        infile >> s1;
        CHECK(s1.compare(">2") == 0);
        infile >> s1;
        CHECK(s1.compare("CTATAGTG") == 0);
        CHECK(std::filesystem::remove("test-write-fasta.fasta"));
    }
    SUBCASE("fst") {
        coati::data_t fasta("", {"1", "2"});
        fasta.out_file.path = "test-write-fasta-fst.fasta";

        VectorFstStdArc fst_write;
        fst_write.AddState();
        fst_write.SetStart(0);
        add_arc(fst_write, 0, 1, 2, 2);  // C -> C
        add_arc(fst_write, 1, 2, 4, 4);  // T -> T
        add_arc(fst_write, 2, 3, 0, 2);  // - -> C
        add_arc(fst_write, 3, 4, 1, 0);  // A -> -
        fst_write.SetFinal(4, 0.0);

        REQUIRE(write_fasta(fasta, fst_write));

        std::ifstream infile("test-write-fasta-fst.fasta");
        std::string s1;
        infile >> s1;
        CHECK(s1.compare(">1") == 0);
        infile >> s1;
        CHECK(s1.compare("CT-A") == 0);
        infile >> s1;
        CHECK(s1.compare(">2") == 0);
        infile >> s1;
        CHECK(s1.compare("CTC-") == 0);
        CHECK(std::filesystem::remove("test-write-fasta-fst.fasta"));
    }
    SUBCASE("data_t size fails - diff number of names and seqs") {
        coati::data_t fasta("", {"1", "2", "3"}, {"CTC", "CTA"});
        CHECK_THROWS_AS(write_fasta(fasta), std::invalid_argument);
    }
}
// GCOVR_EXCL_STOP
}  // namespace coati
