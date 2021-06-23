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

#ifndef DNA_SYMS_HPP
#define DNA_SYMS_HPP

#include <fst/fstlib.h>

inline void fill_symbol_table(fst::SymbolTable& dna_syms) {
    dna_syms.AddSymbol("<eps>", 0);
    dna_syms.AddSymbol("A", 1);
    dna_syms.AddSymbol("C", 2);
    dna_syms.AddSymbol("G", 3);
    dna_syms.AddSymbol("T", 4);
    dna_syms.AddSymbol("U", 4);
    dna_syms.AddSymbol("N", 5);

    dna_syms.AddSymbol("AAA", 11);
    dna_syms.AddSymbol("AAC", 12);
    dna_syms.AddSymbol("AAG", 13);
    dna_syms.AddSymbol("AAT", 14);
    dna_syms.AddSymbol("ACA", 15);
    dna_syms.AddSymbol("ACC", 16);
    dna_syms.AddSymbol("ACG", 17);
    dna_syms.AddSymbol("ACT", 18);
    dna_syms.AddSymbol("AGA", 19);
    dna_syms.AddSymbol("AGC", 20);
    dna_syms.AddSymbol("AGG", 21);
    dna_syms.AddSymbol("AGT", 22);
    dna_syms.AddSymbol("ATA", 23);
    dna_syms.AddSymbol("ATC", 24);
    dna_syms.AddSymbol("ATG", 25);
    dna_syms.AddSymbol("ATT", 26);
    dna_syms.AddSymbol("CAA", 27);
    dna_syms.AddSymbol("CAC", 28);
    dna_syms.AddSymbol("CAG", 29);
    dna_syms.AddSymbol("CAT", 30);
    dna_syms.AddSymbol("CCA", 31);
    dna_syms.AddSymbol("CCC", 32);
    dna_syms.AddSymbol("CCG", 33);
    dna_syms.AddSymbol("CCT", 34);
    dna_syms.AddSymbol("CGA", 35);
    dna_syms.AddSymbol("CGC", 36);
    dna_syms.AddSymbol("CGG", 37);
    dna_syms.AddSymbol("CGT", 38);
    dna_syms.AddSymbol("CTA", 39);
    dna_syms.AddSymbol("CTC", 40);
    dna_syms.AddSymbol("CTG", 41);
    dna_syms.AddSymbol("CTT", 42);
    dna_syms.AddSymbol("GAA", 43);
    dna_syms.AddSymbol("GAC", 44);
    dna_syms.AddSymbol("GAG", 45);
    dna_syms.AddSymbol("GAT", 46);
    dna_syms.AddSymbol("GCA", 47);
    dna_syms.AddSymbol("GCC", 48);
    dna_syms.AddSymbol("GCG", 49);
    dna_syms.AddSymbol("GCT", 50);
    dna_syms.AddSymbol("GGA", 51);
    dna_syms.AddSymbol("GGC", 52);
    dna_syms.AddSymbol("GGG", 53);
    dna_syms.AddSymbol("GGT", 54);
    dna_syms.AddSymbol("GTA", 55);
    dna_syms.AddSymbol("GTC", 56);
    dna_syms.AddSymbol("GTG", 57);
    dna_syms.AddSymbol("GTT", 58);
    dna_syms.AddSymbol("TAA", 59);
    dna_syms.AddSymbol("TAC", 60);
    dna_syms.AddSymbol("TAG", 61);
    dna_syms.AddSymbol("TAT", 62);
    dna_syms.AddSymbol("TCA", 63);
    dna_syms.AddSymbol("TCC", 64);
    dna_syms.AddSymbol("TCG", 65);
    dna_syms.AddSymbol("TCT", 66);
    dna_syms.AddSymbol("TGA", 67);
    dna_syms.AddSymbol("TGC", 68);
    dna_syms.AddSymbol("TGG", 69);
    dna_syms.AddSymbol("TGT", 70);
    dna_syms.AddSymbol("TTA", 71);
    dna_syms.AddSymbol("TTC", 72);
    dna_syms.AddSymbol("TTG", 73);
    dna_syms.AddSymbol("TTT", 74);
}

#endif
