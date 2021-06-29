/*
# Copyright (c) 2020 Reed A. Cartwright <reed@cartwright.ht>
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

#include <charconv>
#include <coati/coati.hpp>
#include <minion.hpp>

#include "verb.hpp"

namespace {
const char *base58_alphabet =
    "123456789ABCDEFGHJKLMNPQRSTUVWXYZabcdefghijkmnopqrstuvwxyz";

std::string base58_encode(uint64_t u) {
    std::string buffer(11, base58_alphabet[0]);
    for(int i = 0; i < 11 && u != 0; ++i) {
        buffer[10 - i] = base58_alphabet[u % 58];
        u = u / 58;
    }
    return buffer;
}

}  // namespace

int main(int argc, char *argv[]) {
    {
        auto check = coati::verb::check_version_number();
        if(check != EXIT_SUCCESS) {
            return check;
        }
    }

    // User-specified Seeds
    std::vector<uint32_t> user_seeds;
    // Go through arguments
    for(int i = 1; i < argc; ++i) {
        int32_t value;
        const char *first = argv[i];
        const char *last = argv[i] + strlen(argv[i]);
        // Does the argument represent a 32-bit signed decimal number
        auto [p, ec] = std::from_chars(first, last, value, 10);
        if(ec == std::errc()) {
            user_seeds.push_back(value);
            continue;
        }
        // Fill 32-bit numbers using the bytes in the string
        for(; first + sizeof(value) < last; first += sizeof(value)) {
            std::memcpy(&value, first, sizeof(value));
            user_seeds.push_back(value);
        }
        // Use the remaining characters
        value = 0;
        std::memcpy(&value, first, last - first);
        user_seeds.push_back(value);
    }

    // Use user-specified seeds or generate a random seed sequence
    auto seeds = (user_seeds.empty())
                     ? minion::create_seed_seq()
                     : minion::SeedSeq32(std::move(user_seeds));

    // Initialize MinionRNG
    minion::Random rand;
    rand.Seed(seeds);

    // Extract State
    minion::Random::state_type state = rand.state();

    // Print
    auto it = state.begin();
    std::cout << base58_encode(*it);
    for(++it; it != state.end(); ++it) {
        std::cout << "-" << base58_encode(*it);
    }
    std::cout << std::endl;

    return EXIT_SUCCESS;
}
