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
#include <random.hpp>

#include "verb.hpp"

int main(int argc, char *argv[]) {
    {
        auto check = coati::verb::check_version_number();
        if(check != EXIT_SUCCESS) {
            return check;
        }
    }

    // Use user-specified seeds or generate a random seed sequence
    auto seeds = (argc < 2)
                     ? fragmites::random::auto_seed_seq()
                     : fragmites::random::string_seed_seq(&argv[1], &argv[argc]);


    // Initialize MinionRNG
    fragmites::random::Random rand;
    rand.Seed(seeds);

    // Extract State
    auto str = fragmites::random::encode_seed(rand.GetSeed());
    // Print State
    std::cout << str << std::endl;

    return EXIT_SUCCESS;
}
