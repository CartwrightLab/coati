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

#include <coati/json.hpp>

namespace coati {
using json_t = nlohmann::json;

/**
 * @brief Convert coati::data_t to json format.
 *
 * @param[in,out] j nlohmann::json json object.
 * @param[in] data coati::data_t object.
 */
void to_json(json_t& j, const data_t& data) {
    for(size_t i = 0; i < data.size(); ++i) {
        j["alignment"][data.names[i]] = data.seqs[i];
    }
    j["score"] = data.score;
}

/**
 * @brief Create coati::data_t object from json input.
 *
 * @param[in,out] data coati::data_t object.
 * @param[in] j nlohmann::json json object.
 */
void from_json(const json_t& j, data_t& data) {
    for(const auto& element : j.at("alignment").items()) {
        data.names.emplace_back(element.key());
        data.seqs.emplace_back(element.value());
    }
    data.score = j.at("score");
}

/**
 * @brief Read json file.
 *
 * @param[in] in std::istream input stream pointing to stdin or file.
 * @param[in] marginal bool true if model is marginal so that FSA are created.
 *
 * @retval coati::data_t content and names of sequences.
 */
coati::data_t read_json(std::istream& in, bool marginal) {
    coati::data_t json;

    // read and convert json format to data_t
    json_t in_json;
    in >> in_json;
    from_json(in_json, json);

    if(!marginal) {  // if model is not marginal, create FSA.
        for(size_t i = 0; i < json.seqs.size(); i++) {
            VectorFstStdArc accept;  // create FSA with sequence
            if(!acceptor(json.seqs[i], accept)) {
                throw std::runtime_error(
                    "Creating acceptor from input json file failed. Exiting!");
            }
            json.fsts.push_back(accept);  // Add FSA
        }
    }

    return json;
}

/// @private
// GCOVR_EXCL_START
TEST_CASE("read_json") {
    std::ofstream outfile;
    std::ifstream in;
    std::string filename{"test-read-json.json"};

    SUBCASE("test-read-json-marg.json") {
        outfile.open(filename);
        REQUIRE(outfile);
        outfile << R"({
  "alignment": {
    "anc": "CTCTGGATAGTC",
    "des": "CTATAGTC"
  },
  "score": 0.1
}
)";
        outfile.close();

        in.open(filename);
        REQUIRE(in);
        coati::data_t data = read_json(in, true);
        REQUIRE(std::filesystem::remove(filename));

        CHECK_EQ(data.names[0], "anc");
        CHECK_EQ(data.names[1], "des");
        CHECK_EQ(data.seqs[0], "CTCTGGATAGTC");
        CHECK_EQ(data.seqs[1], "CTATAGTC");
        CHECK_EQ(data.score, 0.1f);
    }
    SUBCASE("test-read-json-fst.json") {
        outfile.open(filename);
        REQUIRE(outfile);
        outfile << R"({
  "alignment": {
    "anc": "CTCTGGATAGTC",
    "des": "CTATAGTC"
  },
  "score": 0.1
}
)";
        outfile.close();

        in.open(filename);
        REQUIRE(in);
        coati::data_t data = read_json(in, false);
        REQUIRE(std::filesystem::remove(filename));

        CHECK_EQ(data.names[0], "anc");
        CHECK_EQ(data.names[1], "des");
        CHECK_EQ(data.seqs[0], "CTCTGGATAGTC");
        CHECK_EQ(data.seqs[1], "CTATAGTC");
        CHECK_EQ(data.fsts[0].NumStates(), 13);
        CHECK_EQ(data.fsts[1].NumStates(), 9);
        for(int i = 0; i < 12; i++) {
            CHECK_EQ(data.fsts[0].NumArcs(i), 1);
        }

        for(int i = 0; i < 8; i++) {
            CHECK_EQ(data.fsts[0].NumArcs(i), 1);
        }
        CHECK_EQ(data.score, 0.1f);
    }
}
// GCOVR_EXCL_STOP

/**
 * @brief Write content of coati::data_t to json format.
 *
 * @param[in] json coati::data_t object.
 * @param[in] out std::ostream output stream pointing to stdout or file.
 * @param[in] aln coati::VectorFstStdArc FST with alignment path.
 *
 */
void write_json(coati::data_t& json, std::ostream& out,
                const VectorFstStdArc& aln) {
    if(aln.NumStates() > 0) {  // if FST alignment, convert to strings
        coati::utils::fst_to_seqs(json, aln);
    }

    // implicit conversion of data_t to json using `to_json`
    json_t out_json = json;

    out << std::setw(2) << out_json << std::endl;
}

/// @private
// GCOVR_EXCL_START
TEST_CASE("write-json") {
    std::ofstream outfile;
    std::string filename{"test-write-json.json"};
    outfile.open(filename);
    REQUIRE(outfile);
    std::ostream& out = outfile;

    SUBCASE("marginal") {
        coati::data_t json("", {"a", "b"},
                           {"ATGTCTTCTCACAAGACA", "ATGTCTTCTCACAAGACA"});

        write_json(json, out);

        std::ifstream infile(filename);
        std::stringstream ss;
        ss << infile.rdbuf();
        std::string s1 = ss.str();
        CHECK_EQ(s1, R"({
  "alignment": {
    "a": "ATGTCTTCTCACAAGACA",
    "b": "ATGTCTTCTCACAAGACA"
  },
  "score": 0.0
}
)");
        REQUIRE(std::filesystem::remove(filename));
    }
    SUBCASE("FST") {
        coati::data_t json("", {"a", "b"}, {"CT-A", "CTC-"});

        VectorFstStdArc fst_write;
        fst_write.AddState();
        fst_write.SetStart(0);
        add_arc(fst_write, 0, 1, 2, 2);  // C -> C
        add_arc(fst_write, 1, 2, 4, 4);  // T -> T
        add_arc(fst_write, 2, 3, 0, 2);  // - -> C
        add_arc(fst_write, 3, 4, 1, 0);  // A -> -
        fst_write.SetFinal(4, 0.0);

        write_json(json, out, fst_write);

        std::ifstream infile(filename);
        std::stringstream ss;
        ss << infile.rdbuf();
        std::string s1 = ss.str();
        CHECK_EQ(s1, R"({
  "alignment": {
    "a": "CT-A",
    "b": "CTC-"
  },
  "score": 0.0
}
)");
        REQUIRE(std::filesystem::remove(filename));
    }
}
// GCOVR_EXCL_STOP

/**
 * @brief Write content of coati::data_t to json format for coati sample.
 *
 * @param[in] json coati::data_t object.
 * @param[in] out std::ostream output stream pointing to stdout or file.
 * @param[in] iter std::size_t iteration.
 * @param[in] sample_size std::size_t number of iterations.
 *
 */
void write_json(coati::data_t& data, std::ostream& out, size_t iter,
                size_t sample_size) {
    if(iter == 0) {
        out << "[" << std::endl;
    }

    // implicit conversion of data_t to json using `to_json`
    json_t json = data;

    out << std::setw(2) << json;
    if(iter < sample_size - 1) {
        out << "," << std::endl;
    } else {
        out << std::endl << "]" << std::endl;
    }
}

/// @private
// GCOVR_EXCL_START
TEST_CASE("write-json-sample") {
    std::ofstream outfile;
    std::string filename{"test-write-json.json"};
    outfile.open(filename);
    REQUIRE(outfile);
    std::ostream& out = outfile;

    SUBCASE("sample size 1") {
        coati::data_t data("", {"a"}, {"ATGTCTTCTCACAAGACT"});

        write_json(data, out, 0, 1);

        std::ifstream infile(filename);
        std::stringstream ss;
        ss << infile.rdbuf();
        std::string s1 = ss.str();
        CHECK_EQ(s1, R"([
{
  "alignment": {
    "a": "ATGTCTTCTCACAAGACT"
  },
  "score": 0.0
}
]
)");
        REQUIRE(std::filesystem::remove(filename));
    }
    SUBCASE("sample size 2") {
        coati::data_t data("", {"a", "b"},
                           {"ATGTCTTCTCACAAGACT", "ATGTCTTCTCACAAGACA"});

        write_json(data, out, 0, 2);
        write_json(data, out, 1, 2);

        std::ifstream infile(filename);
        std::stringstream ss;
        ss << infile.rdbuf();
        std::string s1 = ss.str();
        CHECK_EQ(s1, R"([
{
  "alignment": {
    "a": "ATGTCTTCTCACAAGACT",
    "b": "ATGTCTTCTCACAAGACA"
  },
  "score": 0.0
},
{
  "alignment": {
    "a": "ATGTCTTCTCACAAGACT",
    "b": "ATGTCTTCTCACAAGACA"
  },
  "score": 0.0
}
]
)");
        REQUIRE(std::filesystem::remove(filename));
    }
}
// GCOVR_EXCL_STOP

}  // namespace coati
