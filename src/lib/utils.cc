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

constexpr size_t PRINT_SIZE = 100;

/* Add arc to FST */
void add_arc(VectorFstStdArc& n2p, int src, int dest, int ilabel, int olabel,
             float weight) {
    if(weight == 1.0) {
        weight = 0.0;
    } else if(weight == 0.0) {
        weight = static_cast<float>(INT_MAX);
    } else {
        weight = -logf(weight);
    }

    if(n2p.NumStates() <= dest) {
        n2p.AddState();
        n2p.AddArc(src, fst::StdArc(ilabel, olabel, weight, dest));
    } else {
        n2p.AddArc(src, fst::StdArc(ilabel, olabel, weight, dest));
    }
}

/* Optimize FST: remove epsilons, determinize, and minimize */
VectorFstStdArc optimize(VectorFstStdArc fst_raw) {
    using fst::StdArc;
    // encode FST
    fst::SymbolTable syms;
    fill_symbol_table(syms);

    fst::EncodeMapper<StdArc> encoder(fst::kEncodeLabels, fst::ENCODE);
    encoder.SetInputSymbols(&syms);
    encoder.SetOutputSymbols(&syms);
    fst::Encode(&fst_raw, &encoder);

    // reduce to more efficient form
    // 1. epsilon removal
    fst::RmEpsilon(&fst_raw);

    // 2. determinize
    VectorFstStdArc fst_det;
    fst::Determinize(fst_raw, &fst_det);

    // 3. minimize
    fst::Minimize(&fst_det);

    // decode
    fst::Decode(&fst_det, encoder);

    return fst_det;
}

/* Read fasta format file */
int read_fasta(fasta_t& fasta_file, std::vector<VectorFstStdArc>& fsts) {
    std::ifstream input(fasta_file.path);
    if(!input.good()) {
        throw std::invalid_argument("Error opening " +
                                    fasta_file.path.string() + ".");
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
                                             fasta_file.path.string() +
                                             " failed. Exiting!");
                }
                fsts.push_back(accept);  // Add FSA
                fasta_file.seq_data.push_back(content);
                name.clear();
            }
            // Add name of sequence
            name = line.substr(1);
            fasta_file.seq_names.push_back(name);
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
            throw std::runtime_error("Creating acceptor from " +
                                     fasta_file.path.string() +
                                     " failed. Exiting!");
        }
        fsts.push_back(accept);
        fasta_file.seq_data.push_back(content);
    }

    return 0;
}

TEST_CASE("read_fasta-fst") {
    // cppcheck-suppress unusedVariable
    std::vector<VectorFstStdArc> fsts;
    std::ofstream outfile;

    SUBCASE("Read test-read-fasta.fasta") {
        outfile.open("test-read-fasta.fasta");
        REQUIRE(outfile);
        outfile << "; comment line" << std::endl;
        outfile << ">1" << std::endl << "CTCTGGATAGTG" << std::endl;
        outfile << ">2" << std::endl << "CTATAGTG" << std::endl;
        outfile.close();

        fasta_t fasta("test-read-fasta.fasta");

        REQUIRE(read_fasta(fasta, fsts) == 0);
        CHECK(std::filesystem::remove(fasta.path));

        CHECK(fasta.seq_data[0] == "CTCTGGATAGTG");
        CHECK(fasta.seq_data[1] == "CTATAGTG");

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
        fasta_t fasta("test-9999999999.fasta");

        REQUIRE_THROWS_AS(read_fasta(fasta, fsts), std::invalid_argument);
    }
}

/* Read fasta format file */
int read_fasta(fasta_t& fasta_file) {
    std::ifstream input(fasta_file.path);
    if(!input.good()) {
        throw std::invalid_argument("Error opening " +
                                    fasta_file.path.string() + ".");
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
                fasta_file.seq_data.push_back(content);
                name.clear();
            }
            // Add name of sequence
            name = line.substr(1);
            fasta_file.seq_names.push_back(name);
            content.clear();
        } else if(!name.empty()) {
            // Remove spaces
            line.erase(remove_if(line.begin(), line.end(), ::isspace),
                       line.end());
            content += line;
        }
    }
    if(!name.empty()) {  // Add last sequence FSA if needed
        fasta_file.seq_data.push_back(content);
    }

    return 0;
}

TEST_CASE("read_fasta") {
    std::ofstream outfile;

    SUBCASE("Read test-read-fasta.fasta") {
        outfile.open("test-read-fasta.fasta");
        REQUIRE(outfile);
        outfile << "; comment line" << std::endl;
        outfile << ">1" << std::endl << "CTCTGGATAGTC" << std::endl;
        outfile << ">2" << std::endl << "CTATAGTC" << std::endl;
        outfile.close();

        fasta_t fasta("test-read-fasta.fasta");

        REQUIRE(read_fasta(fasta) == 0);
        CHECK(std::filesystem::remove(fasta.path));

        CHECK(fasta.seq_names[0] == "1");
        CHECK(fasta.seq_names[1] == "2");
        CHECK(fasta.seq_data[0] == "CTCTGGATAGTC");
        CHECK(fasta.seq_data[1] == "CTATAGTC");
    }

    SUBCASE("Error opening fasta") {
        fasta_t fasta("test-9999999999.fasta");

        REQUIRE_THROWS_AS(read_fasta(fasta), std::invalid_argument);
    }
}

/* Read fasta formatted file with a pair of sequences */
void read_fasta_pair(fasta_t& fasta_file, std::vector<VectorFstStdArc>& fsts,
                     bool fst) {
    if(fst) {
        if(read_fasta(fasta_file, fsts) != 0) {
            throw std::invalid_argument("Error reading " +
                                        fasta_file.path.string() +
                                        " file. Exiting!");
        }
        if(fasta_file.seq_names.size() != 2 ||
           fasta_file.seq_names.size() != fsts.size()) {
            throw std::invalid_argument(
                "Exactly two sequences required. Exiting!");
        }
    } else {
        if(read_fasta(fasta_file) != 0) {
            throw std::invalid_argument("Error reading " +
                                        fasta_file.path.string() +
                                        " file. Exiting!");
        }
        if(fasta_file.seq_names.size() != 2) {
            throw std::invalid_argument(
                "Exactly two sequences required. Exiting!");
        }
    }
}

TEST_CASE("read_fasta_pair") {
    std::ofstream outfile;

    SUBCASE("Read test-no-fst") {
        std::vector<VectorFstStdArc> fsts;
        outfile.open("test-no-fst.fasta");
        REQUIRE(outfile);
        outfile << "; comment line" << std::endl;
        outfile << ">1" << std::endl << "CTCTGGATAGTC" << std::endl;
        outfile << ">2" << std::endl << "CTATAGTC" << std::endl;
        outfile.close();

        fasta_t fasta("test-no-fst.fasta");

        read_fasta_pair(fasta, fsts, false);
        CHECK(std::filesystem::remove(fasta.path));

        CHECK(fasta.seq_names.size() == 2);
    }

    SUBCASE("Read test-fst") {
        std::vector<VectorFstStdArc> fsts;
        outfile.open("test-fst.fasta");
        REQUIRE(outfile);
        outfile << "; comment line" << std::endl;
        outfile << ">1" << std::endl << "CTCTGGATAGTC" << std::endl;
        outfile << ">2" << std::endl << "CTATAGTC" << std::endl;
        outfile.close();

        fasta_t fasta("test-fst.fasta");

        read_fasta_pair(fasta, fsts, true);
        CHECK(std::filesystem::remove(fasta.path));

        CHECK(fasta.seq_names.size() == 2);
        CHECK(fasta.seq_names.size() == fsts.size());
    }

    SUBCASE("Error test-no-fst.fasta") {
        std::vector<VectorFstStdArc> fsts;
        outfile.open("test-no-fst.fasta");
        REQUIRE(outfile);
        outfile << "; comment line" << std::endl;
        outfile << ">1" << std::endl << "CTCTGGATAGTC" << std::endl;
        outfile.close();

        fasta_t fasta("test-no-fst.fasta");
        CHECK(std::filesystem::remove(fasta.path));

        REQUIRE_THROWS_AS(read_fasta_pair(fasta, fsts, false),
                          std::invalid_argument);
    }

    SUBCASE("Error test-fst.fasta") {
        std::vector<VectorFstStdArc> fsts;
        outfile.open("test-fst.fasta");
        REQUIRE(outfile);
        outfile << "; comment line" << std::endl;
        outfile << ">1" << std::endl << "CTCTGGATAGTC" << std::endl;
        outfile << ">2" << std::endl << "CTATAGTC" << std::endl;
        outfile << ">3" << std::endl << "CTATAGTC" << std::endl;
        outfile.close();

        fasta_t fasta("test-fst.fasta");
        CHECK(std::filesystem::remove(fasta.path));

        REQUIRE_THROWS_AS(read_fasta_pair(fasta, fsts, true),
                          std::invalid_argument);
    }
}

/* Write alignment in Fasta format */
int write_fasta(fasta_t& fasta_file) {
    std::ofstream outfile;
    outfile.open(fasta_file.path);
    if(!outfile) {
        throw std::invalid_argument("Opening output file failed.");
    }

    for(size_t i = 0; i < fasta_file.seq_names.size(); i++) {
        outfile << ">" << fasta_file.seq_names[i] << std::endl
                << fasta_file.seq_data[i] << std::endl;
    }
    outfile.close();

    return EXIT_SUCCESS;
}

TEST_CASE("write_fasta") {
    SUBCASE("Write fasta") {
        fasta_t fasta("test-write-fasta.fasta", {"1", "2"},
                      {"CTCTGGATAGTG", "CTATAGTG"});

        REQUIRE(write_fasta(fasta) == 0);

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
int write_fasta(VectorFstStdArc& aln, fasta_t& fasta_file) {
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

    fasta_file.seq_data.push_back(seq1);
    fasta_file.seq_data.push_back(seq2);

    // map all epsilons (<eps>) to gaps (-)
    while(fasta_file.seq_data[0].find("<eps>") != std::string::npos) {
        fasta_file.seq_data[0].replace(fasta_file.seq_data[0].find("<eps>"), 5,
                                       "-");
    }
    while(fasta_file.seq_data[1].find("<eps>") != std::string::npos) {
        fasta_file.seq_data[1].replace(fasta_file.seq_data[1].find("<eps>"), 5,
                                       "-");
    }

    // write aligned sequences to file
    write_fasta(fasta_file);
    return EXIT_SUCCESS;
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

    REQUIRE(write_fasta(fst_write, fasta) == 0);

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

/* Write alignment in PHYLIP format */
int write_phylip(fasta_t& fasta_file) {
    std::ofstream outfile;
    outfile.open(fasta_file.path);
    if(!outfile) {
        throw std::invalid_argument("Opening output file failed.");
    }

    // write aligned sequences to file
    outfile << fasta_file.seq_names.size() << " "
            << fasta_file.seq_data[0].length() << std::endl;
    size_t i = PRINT_SIZE - 4 -
               std::max(fasta_file.seq_names[0].length(),
                        fasta_file.seq_names[1].length());
    for(size_t j = 0; j < fasta_file.seq_names.size(); j++) {
        outfile << fasta_file.seq_names[j] << "\t"
                << fasta_file.seq_data[j].substr(0, i) << std::endl;
    }
    outfile << std::endl;

    for(; i < fasta_file.seq_data[0].length(); i += PRINT_SIZE) {
        for(size_t j = 0; j < fasta_file.seq_names.size(); j++) {
            outfile << fasta_file.seq_data[j].substr(i, PRINT_SIZE)
                    << std::endl;
        }
        outfile << std::endl;
    }

    return EXIT_SUCCESS;
}

TEST_CASE("write_phylip") {
    SUBCASE("Short sequences") {
        fasta_t fasta("test-write-phylip.phylip", {"1", "2"},
                      {"CTCTGGATAGTG", "CT----ATAGTG"});

        REQUIRE(write_phylip(fasta) == 0);

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
        fasta_t fasta("test-write-phylip.phylip", {"1", "2"},
                      {"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
                       "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
                       "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
                       "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"});

        REQUIRE(write_phylip(fasta) == 0);

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
        fasta_t fasta("", {"1", "2"}, {"CTCTGGATAGTG", "CTATAGTG"});

        REQUIRE_THROWS_AS(write_phylip(fasta), std::invalid_argument);
    }
}

/* Write shortest path (alignment) in PHYLIP format */
int write_phylip(VectorFstStdArc& aln, fasta_t& fasta_file) {
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

    fasta_file.seq_data.push_back(seq1);
    fasta_file.seq_data.push_back(seq2);

    // map all epsilons (<eps>) to gaps (-)
    while(fasta_file.seq_data[0].find("<eps>") != std::string::npos) {
        fasta_file.seq_data[0].replace(fasta_file.seq_data[0].find("<eps>"), 5,
                                       "-");
    }
    while(fasta_file.seq_data[1].find("<eps>") != std::string::npos) {
        fasta_file.seq_data[1].replace(fasta_file.seq_data[1].find("<eps>"), 5,
                                       "-");
    }

    return write_phylip(fasta_file);
}

TEST_CASE("write_phylip-fst") {
    fasta_t fasta("test-write-phylip.phylip", {"1", "2"});

    VectorFstStdArc fst_write;
    fst_write.AddState();
    fst_write.SetStart(0);
    add_arc(fst_write, 0, 1, 2, 2);  // C -> C
    add_arc(fst_write, 1, 2, 4, 4);  // T -> T
    add_arc(fst_write, 2, 3, 0, 2);  // - -> C
    add_arc(fst_write, 3, 4, 1, 0);  // A -> -
    fst_write.SetFinal(4, 0.0);

    REQUIRE(write_phylip(fst_write, fasta) == 0);
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

/* Create FSAs (acceptors) from a fasta file*/
bool acceptor(std::string content, VectorFstStdArc& accept) {
    std::map<char, int> syms = {{'-', 0}, {'A', 1}, {'C', 2}, {'G', 3},
                                {'T', 4}, {'U', 4}, {'N', 5}};

    // Add initial state
    accept.AddState();
    accept.SetStart(0);

    for(int i = 0, content_len = static_cast<int>(content.length());
        i < content_len; i++) {
        add_arc(accept, i, i + 1, syms.at(content[i]), syms.at(content[i]));
    }

    // Add final state and run an FST sanity check (Verify)
    accept.SetFinal(accept.NumStates() - 1, 0.0);
    return Verify(accept);
}

/* Hamming distance between two codons */
int cod_distance(uint8_t cod1, uint8_t cod2) {
    int distance = 0;

    distance += (((cod1 & 48) >> 4) == ((cod2 & 48) >> 4) ? 0 : 1);
    distance += (((cod1 & 12) >> 2) == ((cod2 & 12) >> 2) ? 0 : 1);
    distance += ((cod1 & 3) == (cod2 & 3) ? 0 : 1);

    return distance;
}

/* Cast codon to position in codon list AAA->0, AAAC->1 ... TTT->63 */
int cod_int(std::string codon) {
    unsigned char pos0 = codon[0];
    unsigned char pos1 = codon[1];
    unsigned char pos2 = codon[2];

    return (nt4_table[pos0] << 4) | (nt4_table[pos1] << 2) | nt4_table[pos2];
}

/* Read substitution rate matrix from a CSV file */
Matrix parse_matrix_csv(const std::string& file) {
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

    Matrix P(64, 64, Q);

    return P;
}

TEST_CASE("parse_matrix_csv") {
    std::ofstream outfile;
    Matrix P(mg94_p(0.0133, 0.2));

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
    Matrix P_test(parse_matrix_csv("test-marg-matrix.csv"));
    for(auto i = 0; i < 64; i++) {
        for(auto j = 0; j < 64; j++) {
            CHECK(P(i, j) == doctest::Approx(P_test(i, j)));
        }
    }
    CHECK(std::filesystem::remove("test-marg-matrix.csv"));
}
