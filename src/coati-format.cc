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

#include <CLI11.hpp>
#include <coati/format.hpp>
#include <coati/io.hpp>
#include <coati/utils.hpp>

int main(int argc, char* argv[]) {
    coati::args_t args;

    // Parse command line options
    CLI::App format{
        "coati format - convert between formats, extract and/or reoder "
        "sequences\n"};
    coati::utils::set_options_format(format, args);
    CLI11_PARSE(format, argc, argv);

    // if no input specified, use cin and json format as default
    if(args.aln.data.path.empty()) {
        args.aln.data.path = "json:-";
    }

    try {
        // read input data
        args.aln.data = coati::io::read_input(args.aln);

        return coati::format_sequences(args.format, args.aln);
    } catch(const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << std::endl;
    }
    return EXIT_FAILURE;
}
