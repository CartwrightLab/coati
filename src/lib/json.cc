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
    j["data"]["names"] = data.names;
    j["data"]["seqs"] = data.seqs;
}

/**
 * @brief Create coati::data_t object from json input.
 *
 * @param[in,out] data coati::data_t object.
 * @param[in] j nlohmann::json json object.
 */
void from_json(const json_t& j, data_t& data) {
    for(size_t i = 0; i < j.at("data").size(); i++) {
        data.names.emplace_back(j.at("data").at("names")[i]);
        data.seqs.emplace_back(j.at("data").at("seqs")[i]);
    }
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
        outfile << "{\"data\":{\"names\" : [\"anc\",\"des\"], \"seqs\" : "
                   "[\"CTCTGGATAGTC\",\"CTATAGTC\"]}}"
                << std::endl;
        outfile.close();

        in.open(filename);
        REQUIRE(in);
        coati::data_t json = read_json(in, true);
        REQUIRE(std::filesystem::remove(filename));

        CHECK_EQ(json.names[0], "anc");
        CHECK_EQ(json.names[1], "des");
        CHECK_EQ(json.seqs[0], "CTCTGGATAGTC");
        CHECK_EQ(json.seqs[1], "CTATAGTC");
    }
    SUBCASE("test-read-json-fst.json") {
        outfile.open(filename);
        REQUIRE(outfile);
        outfile << "{\"data\":{\"names\" : [\"anc\",\"des\"], \"seqs\" : "
                   "[\"CTCTGGATAGTC\",\"CTATAGTC\"]}}"
                << std::endl;
        outfile.close();

        in.open(filename);
        REQUIRE(in);
        coati::data_t json = read_json(in, false);
        REQUIRE(std::filesystem::remove(filename));

        CHECK_EQ(json.names[0], "anc");
        CHECK_EQ(json.names[1], "des");
        CHECK_EQ(json.seqs[0], "CTCTGGATAGTC");
        CHECK_EQ(json.seqs[1], "CTATAGTC");
        CHECK_EQ(json.fsts[0].NumStates(), 13);
        CHECK_EQ(json.fsts[1].NumStates(), 9);
        for(int i = 0; i < 12; i++) {
            CHECK_EQ(json.fsts[0].NumArcs(i), 1);
        }

        for(int i = 0; i < 8; i++) {
            CHECK_EQ(json.fsts[0].NumArcs(i), 1);
        }
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
void write_json(coati::data_t& json, std::ostream& out) {
    // implicit conversion of data_t to json using `to_json`
    json_t out_json = json;

    out << out_json << std::endl;
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
                           {"ATGTCTTCTCACAAGACT", "ATGTCTTCTCACAAGACT"});

        write_json(json, out);

        std::ifstream infile(filename);
        std::string s1;
        infile >> s1;
        CHECK_EQ(s1,
                 "{\"data\":{\"names\":[\"a\",\"b\"],\"seqs\":"
                 "[\"ATGTCTTCTCACAAGACT\",\"ATGTCTTCTCACAAGACT\"]}}");
        REQUIRE(std::filesystem::remove(filename));
    }
}
// GCOVR_EXCL_STOP

}  // namespace coati
