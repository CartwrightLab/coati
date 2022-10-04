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
 * \brief Convert coati::data_t to json format.
 *
 * @param[in,out] j nlohmann::json json object.
 * @param[in] data coati::data_t object.
 */
void to_json(json_t& j, const data_t& data) {
    j["data"]["names"] = data.names;
    j["data"]["seqs"] = data.seqs;
}

/**
 * \brief Create coati::data_t object from json input.
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
 * \brief Read json file.
 *
 * @param[in] f_path std::string path to input file.
 * @param[in] marginal bool true if model is marginal so that FSA are created.
 *
 * \return coati::data_t object.
 */
coati::data_t read_json(const std::string& f_path, bool marginal) {
    coati::data_t json(f_path);

    // set input pointer and file type
    coati::file_type_t in_type = coati::utils::extract_file_type(f_path);
    std::istream* is = coati::io::set_istream(in_type.path);
    std::istream& in = *is;

    // read and convert json format to data_t
    json_t in_json;
    in >> in_json;
    from_json(in_json, json);

    if(!marginal) {  // if model is not marginal, create FSA.
        for(size_t i = 0; i < json.seqs.size(); i++) {
            VectorFstStdArc accept;  // create FSA with sequence
            if(!acceptor(json.seqs[i], accept)) {
                throw std::runtime_error("Creating acceptor from " + f_path +
                                         " failed. Exiting!");
            }
            json.fsts.push_back(accept);  // Add FSA
        }
    }

    return json;
}

/// @private
// GCOVR_EXCL_START
TEST_CASE("read-json") {
    std::ofstream outfile;

    SUBCASE("test-read-json-marg.json") {
        std::string filename{"test-read-json-marg.json"};
        outfile.open(filename);
        REQUIRE(outfile);
        outfile << "{\"data\":{\"names\" : [\"anc\",\"des\"], \"seqs\" : "
                   "[\"CTCTGGATAGTC\",\"CTATAGTC\"]}}"
                << std::endl;
        outfile.close();

        coati::data_t json = read_json(filename, true);
        CHECK(std::filesystem::remove(json.path));

        CHECK(json.names[0] == "anc");
        CHECK(json.names[1] == "des");
        CHECK(json.seqs[0] == "CTCTGGATAGTC");
        CHECK(json.seqs[1] == "CTATAGTC");
    }
    SUBCASE("test-read-json-fst.json") {
        std::string filename{"test-read-json-fst.json"};
        outfile.open(filename);
        REQUIRE(outfile);
        outfile << "{\"data\":{\"names\" : [\"anc\",\"des\"], \"seqs\" : "
                   "[\"CTCTGGATAGTC\",\"CTATAGTC\"]}}"
                << std::endl;
        outfile.close();

        coati::data_t json = read_json(filename, false);
        CHECK(std::filesystem::remove(json.path));

        CHECK(json.names[0] == "anc");
        CHECK(json.names[1] == "des");
        CHECK(json.seqs[0] == "CTCTGGATAGTC");
        CHECK(json.seqs[1] == "CTATAGTC");
        CHECK(json.fsts[0].NumStates() == 13);
        CHECK(json.fsts[1].NumStates() == 9);
        for(int i = 0; i < 12; i++) {
            CHECK(json.fsts[0].NumArcs(i) == 1);
        }

        for(int i = 0; i < 8; i++) {
            CHECK(json.fsts[0].NumArcs(i) == 1);
        }
    }
    SUBCASE("invalid input argument") {
        std::string filename{"test-read-json.json"};
        CHECK_THROWS_AS(coati::read_json(filename, true),
                        std::invalid_argument);
    }
}
// GCOVR_EXCL_STOP

/**
 * \brief Write content of coati::data_t to json format.
 *
 * @param[in] json coati::data_t object.
 * @param[in] aln coati::VectorFstStdArc FST with alignment path.
 *
 */
void write_json(coati::data_t& json, const VectorFstStdArc& aln) {
    if(aln.NumStates() > 0) {  // if FST alignment
        coati::utils::fst_to_seqs(json, aln);
    }
    // set output pointer
    std::ostream* os = coati::io::set_ostream(json.out_file.path);
    std::ostream& out = *os;

    json_t out_json = json;

    out << out_json << std::endl;
}

/// @private
// GCOVR_EXCL_START
TEST_CASE("write-json") {
    SUBCASE("test-write-json-marg.json") {
        coati::data_t json("", {"a", "b"},
                           {"ATGTCTTCTCACAAGACT", "ATGTCTTCTCACAAGACT"});
        json.out_file.path = "test-write-json-marg.json";

        write_json(json);

        std::ifstream infile("test-write-json-marg.json");
        std::string s1;
        infile >> s1;
        CHECK(s1 ==
              "{\"data\":{\"names\":[\"a\",\"b\"],\"seqs\":"
              "[\"ATGTCTTCTCACAAGACT\",\"ATGTCTTCTCACAAGACT\"]}}");
        CHECK(std::filesystem::remove(json.out_file.path));
    }
    SUBCASE("test-write-json-fst.json") {
        coati::data_t json("", {"a", "b"}, {"CT-A", "CTC-"});
        json.out_file.path = "test-write-json-fst.json";

        VectorFstStdArc fst_write;
        fst_write.AddState();
        fst_write.SetStart(0);
        add_arc(fst_write, 0, 1, 2, 2);  // C -> C
        add_arc(fst_write, 1, 2, 4, 4);  // T -> T
        add_arc(fst_write, 2, 3, 0, 2);  // - -> C
        add_arc(fst_write, 3, 4, 1, 0);  // A -> -
        fst_write.SetFinal(4, 0.0);

        write_json(json, fst_write);

        std::ifstream infile(json.out_file.path);
        std::string s1;
        infile >> s1;
        CHECK(s1 ==
              "{\"data\":{\"names\":[\"a\",\"b\"],\"seqs\":"
              "[\"CT-A\",\"CTC-\"]}}");
        CHECK(std::filesystem::remove(json.out_file.path));
    }
}
// GCOVR_EXCL_STOP

}  // namespace coati
