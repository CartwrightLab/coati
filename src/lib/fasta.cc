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

/* Read fasta format file */
fasta_t read_fasta(const std::string& f_path,
                   std::vector<VectorFstStdArc>& fsts) {
    fasta_t fasta(f_path);
    std::ifstream input(f_path);
    if(!input.good()) {
        throw std::invalid_argument("Error opening " + f_path + ".");
    }

    std::string line, name, content;
    while(getline(input, line).good()) {
        if(line.empty()) {
            continue;  // omit empty lines
        }
        if(line[0] == ';') {
            continue;
        }
        if(line[0] == '>') {  // Identifier marker
            if(!name.empty()) {
                VectorFstStdArc accept;  // create FSA with sequence
                if(!acceptor(content, accept)) {
                    throw std::runtime_error("Creating acceptor from " +
                                             f_path + " failed. Exiting!");
                }
                fsts.push_back(accept);  // Add FSA
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
        VectorFstStdArc accept;
        if(!acceptor(content, accept)) {
            throw std::runtime_error("Creating acceptor from " + f_path +
                                     " failed. Exiting!");
        }
        fsts.push_back(accept);
        fasta.seqs.push_back(content);
    }

    return fasta;
}

TEST_CASE("read_fasta-fst") {
    // cppcheck-suppress unusedVariable
    std::vector<VectorFstStdArc> fsts;
    std::ofstream outfile;

    SUBCASE("Read test-read-fasta.fasta") {
        std::string filename{"test-read-fasta.fasta"};
        outfile.open(filename);
        REQUIRE(outfile);
        outfile << "; comment line" << std::endl;
        outfile << ">1" << std::endl << "CTCTGGATAGTG" << std::endl;
        outfile << ">2" << std::endl << "CTATAGTG" << std::endl;
        outfile.close();

        fasta_t fasta = read_fasta(filename, fsts);
        CHECK(std::filesystem::remove(fasta.path.string()));

        CHECK(fasta.seqs[0] == "CTCTGGATAGTG");
        CHECK(fasta.seqs[1] == "CTATAGTG");

        CHECK(fsts[0].NumStates() == 13);
        CHECK(fsts[1].NumStates() == 9);

        for(int i = 0; i < 12; i++) {
            CHECK(fsts[0].NumArcs(i) == 1);
        }

        for(int i = 0; i < 8; i++) {
            CHECK(fsts[0].NumArcs(i) == 1);
        }
    }

    SUBCASE("Error opening fasta") {
        REQUIRE_THROWS_AS(read_fasta("test-9999999999.fasta", fsts),
                          std::invalid_argument);
    }
}

/* Read fasta format file */
fasta_t read_fasta(const std::string& f_path) {
    fasta_t fasta(f_path);
    std::ifstream input(f_path);
    if(!input.good()) {
        throw std::invalid_argument("Error opening " + f_path + ".");
    }

    std::string line, name, content;
    while(getline(input, line).good()) {
        if(line.empty()) {
            continue;  // omit empty lines
        }
        if(line[0] == ';') {
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

    return fasta;
}

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

        fasta_t fasta = read_fasta(filename);
        CHECK(std::filesystem::remove(fasta.path));

        CHECK(fasta.names[0] == "1");
        CHECK(fasta.names[1] == "2");
        CHECK(fasta.seqs[0] == "CTCTGGATAGTC");
        CHECK(fasta.seqs[1] == "CTATAGTC");
    }

    SUBCASE("Error opening fasta") {
        REQUIRE_THROWS_AS(read_fasta("test-9999999999.fasta"),
                          std::invalid_argument);
    }
}

/* Write alignment in Fasta format */
bool write_fasta(const fasta_t& fasta) {
    std::ofstream outfile;
    outfile.open(fasta.path);
    if(!outfile) {
        throw std::invalid_argument("Opening output file" +
                                    fasta.path.string() + "  failed.");
    }

    for(size_t i = 0; i < fasta.size(); i++) {
        outfile << ">" << fasta.names[i] << std::endl
                << fasta.seqs[i] << std::endl;
    }
    outfile.close();

    return true;
}

TEST_CASE("write_fasta") {
    SUBCASE("Write fasta") {
        fasta_t fasta("test-write-fasta.fasta", {"1", "2"},
                      {"CTCTGGATAGTG", "CTATAGTG"});

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

    SUBCASE("Opening file fails") {
        fasta_t f("", {"1", "2"}, {"CTCTGGATAGTG", "CTATAGTG"});

        REQUIRE_THROWS_AS(write_fasta(f), std::invalid_argument);
    }
}

/* Write shortest path (alignment) in Fasta format */
bool write_fasta(const VectorFstStdArc& aln, fasta_t& fasta) {
    fst::SymbolTable symbols;
    fill_symbol_table(symbols);

    std::string seq1, seq2;
    fst::StateIterator<fst::StdFst> siter(aln);  // FST state iterator
    for(int i = 0; i < (aln.NumStates() - 1); siter.Next(), i++) {
        fst::ArcIteratorData<fst::StdArc> data;
        aln.InitArcIterator(siter.Value(), &data);
        seq1.append(symbols.Find(data.arcs[0].ilabel));
        seq2.append(symbols.Find(data.arcs[0].olabel));
    }

    fasta.seqs.resize(2);
    fasta.seqs[0] = seq1;
    fasta.seqs[1] = seq2;

    // map all epsilons (<eps>) to gaps (-)
    while(fasta.seqs[0].find("<eps>") != std::string::npos) {
        fasta.seqs[0].replace(fasta.seqs[0].find("<eps>"), 5, "-");
    }
    while(fasta.seqs[1].find("<eps>") != std::string::npos) {
        fasta.seqs[1].replace(fasta.seqs[1].find("<eps>"), 5, "-");
    }

    // write aligned sequences to file
    write_fasta(fasta);
    return true;
}

TEST_CASE("write_fasta-fst") {
    fasta_t fasta("test-write-fasta.fasta", {"1", "2"});

    VectorFstStdArc fst_write;
    fst_write.AddState();
    fst_write.SetStart(0);
    add_arc(fst_write, 0, 1, 2, 2);  // C -> C
    add_arc(fst_write, 1, 2, 4, 4);  // T -> T
    add_arc(fst_write, 2, 3, 0, 2);  // - -> C
    add_arc(fst_write, 3, 4, 1, 0);  // A -> -
    fst_write.SetFinal(4, 0.0);

    REQUIRE(write_fasta(fst_write, fasta));

    std::ifstream infile("test-write-fasta.fasta");
    std::string s1;
    infile >> s1;
    CHECK(s1.compare(">1") == 0);
    infile >> s1;
    CHECK(s1.compare("CT-A") == 0);
    infile >> s1;
    CHECK(s1.compare(">2") == 0);
    infile >> s1;
    CHECK(s1.compare("CTC-") == 0);
    CHECK(std::filesystem::remove("test-write-fasta.fasta"));
}
}  // namespace coati
