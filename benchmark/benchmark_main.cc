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

#include <benchmark/benchmark.h>

#include <coati/align_marginal.hpp>
#include <fasta.tcc>
#include <string>

static void BM_marg_alignment_short(benchmark::State& state) {
    coati::utils::args_t args("", {"1", "2"}, {"CTCTGGATAGTG", "CTATAGTG"},
                              "m-coati", "", "benchmark-marg-alignment.fasta");
    coati::utils::alignment_t aln;
    aln.fasta.path = args.output;
    aln.fasta.names = args.fasta.names;
    coati::utils::set_subst(args, aln);

    for(auto _ : state) {
        coati::marg_alignment(args, aln);
    }

    std::filesystem::remove(aln.fasta.path);
}

BENCHMARK(BM_marg_alignment_short);

static void BM_marg_alignment_1k(benchmark::State& state) {
    std::string seq_a(example_1k_a, 1002);
    std::string seq_b(example_1k_b, 1002);
    coati::utils::args_t args("", {"A", "B"}, {seq_a, seq_b}, "m-coati", "",
                              "benchmark-marg-alignment.fasta");
    coati::utils::alignment_t aln;
    aln.fasta.path = args.output;
    aln.fasta.names = args.fasta.names;
    coati::utils::set_subst(args, aln);

    for(auto _ : state) {
        coati::marg_alignment(args, aln);
    }

    std::filesystem::remove(aln.fasta.path);
}

BENCHMARK(BM_marg_alignment_1k);

BENCHMARK_MAIN();
