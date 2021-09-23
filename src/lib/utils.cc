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

/* Write alignment in PHYLIP format */
bool write_phylip(coati::fasta_t& fasta) {
    std::ofstream outfile;
    outfile.open(fasta.path);
    if(!outfile) {
        throw std::invalid_argument("Opening output file failed.");
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

/* Write shortest path (alignment) in PHYLIP format */
bool write_phylip(VectorFstStdArc& aln, coati::fasta_t& fasta) {
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

    fasta.seqs.push_back(seq1);
    fasta.seqs.push_back(seq2);

    // map all epsilons (<eps>) to gaps (-)
    while(fasta.seqs[0].find("<eps>") != std::string::npos) {
        fasta.seqs[0].replace(fasta.seqs[0].find("<eps>"), 5, "-");
    }
    while(fasta.seqs[1].find("<eps>") != std::string::npos) {
        fasta.seqs[1].replace(fasta.seqs[1].find("<eps>"), 5, "-");
    }

    return write_phylip(fasta);
}

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

/* Hamming distance between two codons */
int cod_distance(uint8_t cod1, uint8_t cod2) {
    int distance = 0;

    distance += (((cod1 & 48) >> 4) == ((cod2 & 48) >> 4) ? 0 : 1);
    distance += (((cod1 & 12) >> 2) == ((cod2 & 12) >> 2) ? 0 : 1);
    distance += ((cod1 & 3) == (cod2 & 3) ? 0 : 1);

    return distance;
}

/* Cast codon to position in codon list AAA->0, AAAC->1 ... TTT->63 */
int cod_int(const std::string& codon) {
    unsigned char pos0 = codon[0];
    unsigned char pos1 = codon[1];
    unsigned char pos2 = codon[2];

    return (nt4_table[pos0] << 4) | (nt4_table[pos1] << 2) | nt4_table[pos2];
}

/* Read substitution rate matrix from a CSV file */
coati::Matrix<coati::float_t> parse_matrix_csv(const std::string& file) {
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

/* Setup command line options */
void set_cli_options(CLI::App& app, coati::utils::args_t& in_data,
                     const std::string& command) {
    // Add new options/flags
    app.add_option("fasta", in_data.fasta.path, "Fasta file path")
        ->required()
        ->check(CLI::ExistingFile);
    if(command.compare("msa") == 0) {
        app.add_option("tree", in_data.tree, "Newick phylogenetic tree")
            ->required()
            ->check(CLI::ExistingFile);
        app.add_option("reference", in_data.ref, "Reference sequence")
            ->required();
    }
    app.add_option("-m,--model", in_data.model, "Substitution model");
    if(command.compare("alignpair") == 0) {
        app.add_option("-t,--time", in_data.br_len,
                       "Evolutionary time/branch length")
            ->check(CLI::PositiveNumber);
        app.add_option("-l,--weight", in_data.weight_file,
                       "Write alignment score to file");
        app.add_flag("-s,--score", in_data.score, "Score alignment");
    }
    app.add_option("-o,--output", in_data.output, "Alignment output file");
    app.add_option("-g,--gap-open", in_data.gap.open, "Gap opening score")
        ->check(CLI::PositiveNumber);
    app.add_option("-e,--gap-extend", in_data.gap.extend, "Gap extension score")
        ->check(CLI::PositiveNumber);
    app.add_option("-w,--omega", in_data.omega, "Nonsynonymous-synonymous bias")
        ->check(CLI::PositiveNumber);
    app.add_option("-p,--pi", in_data.pi, "Nucleotide frequencies (A C G T)")
        ->expected(4);
    app.add_option("-n,--gap-len", in_data.gap.len, "Set gap unit size");
}

/* Encode ( as vector<unsigned char>) ancestor (ref) sequence as codon & phase,
 *      descendant as nucs */
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
