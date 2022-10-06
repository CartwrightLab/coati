/*
# Copyright (c) 2020-2022 Juan J. Garcia Mesa <juanjosegarciamesa@gmail.com>
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

#include <CLI11.hpp>
#include <coati/align_fst.hpp>
#include <coati/align_marginal.hpp>
#include <coati/io.hpp>
#include <coati/utils.hpp>

int main(int argc, char* argv[]) {
    coati::args_t args;

    // Parse command line options
    CLI::App alignpair;
    coati::utils::set_options_sample(alignpair, args);
    CLI11_PARSE(alignpair, argc, argv);

    if(!args.aln.is_marginal()) {
        std::cerr
            << "ERROR: Sampling only available with models m-coati or m-ecm."
            << std::endl;
        return EXIT_FAILURE;
    }

    // Use user-specified seeds or generate a random seed sequence
    coati::random_t rand;
    auto seeds = args.sample.seeds.empty()
                     ? fragmites::random::auto_seed_seq()
                     : fragmites::random::string_seed_seq(
                           args.sample.seeds.begin(), args.sample.seeds.end());
    rand.Seed(seeds);

    try {
        coati::marg_sample(args.aln, args.sample.sample_size, rand);
    } catch(const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
