/*
# Copyright (c) 2020 Juan J. Garcia Mesa <juanjosegarciamesa@gmail.com>
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

#ifndef ALIGN_H
#define ALIGN_H

// #include <fstream>
#include <coati/mut_models.hpp>
#include <boost/filesystem.hpp>

using namespace fst;
using namespace std;

int mcoati(string fasta, vector<string> seq_names,
	vector<string> sequences, bool score, string weight_f, string output);
int fst_alignment(string model, vector<VectorFst<StdArc>>& fsts, vector<string> seq_names,
	string fasta, string weight_f, string output, vector<string> sequences);

#endif
