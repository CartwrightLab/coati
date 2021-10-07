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
 * \brief Write alignment in PHYLIP format
 *
 * @param[in] fasta coati::fasta_t path, names, and sequences.
 */
bool write_phylip(const coati::fasta_t& fasta) {
    std::ofstream outfile;
    outfile.open(fasta.path);
    if(!outfile) {
        throw std::invalid_argument("Opening output file " +
                                    fasta.path.string() + " failed.");
    }

    // write aligned sequences to file
    outfile << fasta.size() << " " << fasta.seqs[0].length() << std::endl;
    size_t i = PRINT_SIZE - 4 -
               std::max(fasta.names[0].length(), fasta.names[1].length());
    for(size_t j = 0; j < fasta.size(); j++) {
        outfile << fasta.names[j] << "\t" << fasta.seqs[j].substr(0, i)
                << std::endl;
    }
    outfile << std::endl;

    for(; i < fasta.seqs[0].length(); i += PRINT_SIZE) {
        for(size_t j = 0; j < fasta.size(); j++) {
            outfile << fasta.seqs[j].substr(i, PRINT_SIZE) << std::endl;
        }
        outfile << std::endl;
    }

    return true;
}

/// @private
TEST_CASE("write_phylip") {
    SUBCASE("Short sequences") {
        coati::fasta_t fasta("test-write-phylip.phylip", {"1", "2"},
                             {"CTCTGGATAGTG", "CT----ATAGTG"});

        REQUIRE(write_phylip(fasta));

        std::ifstream infile("test-write-phylip.phylip");
        std::string s1, s2;

        infile >> s1 >> s2;
        CHECK(s1.compare("2") == 0);
        CHECK(s2.compare("12") == 0);

        infile >> s1 >> s2;
        CHECK(s1.compare("1") == 0);
        CHECK(s2.compare("CTCTGGATAGTG") == 0);

        infile >> s1 >> s2;
        CHECK(s1.compare("2") == 0);
        CHECK(s2.compare("CT----ATAGTG") == 0);

        CHECK(std::filesystem::remove("test-write-phylip.phylip"));
    }

    SUBCASE("Multi-line sequences") {
        coati::fasta_t fasta(
            "test-write-phylip.phylip", {"1", "2"},
            {"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
             "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
             "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
             "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"});

        REQUIRE(write_phylip(fasta));

        std::ifstream infile("test-write-phylip.phylip");
        std::string s1, s2;

        infile >> s1 >> s2;
        CHECK(s1.compare("2") == 0);
        CHECK(s2.compare("100") == 0);

        infile >> s1 >> s2;
        CHECK(s1.compare("1") == 0);
        CHECK(s2.compare("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
                         "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA") == 0);

        infile >> s1 >> s2;
        CHECK(s1.compare("2") == 0);
        CHECK(s2.compare("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
                         "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA") == 0);

        CHECK(std::filesystem::remove("test-write-phylip.phylip"));
    }

    SUBCASE("Opening file fails") {
        coati::fasta_t fasta("", {"1", "2"}, {"CTCTGGATAGTG", "CTATAGTG"});

        REQUIRE_THROWS_AS(write_phylip(fasta), std::invalid_argument);
    }
}

/**
 * \brief Write shortest path (alignment) in PHYLIP format
 *
 * @param[in] aln coati::VectorFstStdArc FST containing alignment.
 * @param[in] fasta coati::fasta_t path and names.
 *
 */
bool write_phylip(const VectorFstStdArc& aln, coati::fasta_t& fasta) {
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

    return write_phylip(fasta);
}

/// @private
TEST_CASE("write_phylip-fst") {
    coati::fasta_t fasta("test-write-phylip.phylip", {"1", "2"});

    VectorFstStdArc fst_write;
    fst_write.AddState();
    fst_write.SetStart(0);
    add_arc(fst_write, 0, 1, 2, 2);  // C -> C
    add_arc(fst_write, 1, 2, 4, 4);  // T -> T
    add_arc(fst_write, 2, 3, 0, 2);  // - -> C
    add_arc(fst_write, 3, 4, 1, 0);  // A -> -
    fst_write.SetFinal(4, 0.0);

    REQUIRE(write_phylip(fst_write, fasta));
    std::ifstream infile("test-write-phylip.phylip");
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

    CHECK(std::filesystem::remove("test-write-phylip.phylip"));
}
}  // namespace coati
