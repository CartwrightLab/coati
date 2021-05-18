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

#ifndef DNA_SYMS_H
#define DNA_SYMS_H

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

    dna_syms.AddSymbol("AAA1", 101);
    dna_syms.AddSymbol("AAA2", 102);
    dna_syms.AddSymbol("AAA3", 103);
    dna_syms.AddSymbol("AAC1", 104);
    dna_syms.AddSymbol("AAC2", 105);
    dna_syms.AddSymbol("AAC3", 106);
    dna_syms.AddSymbol("AAG1", 107);
    dna_syms.AddSymbol("AAG2", 108);
    dna_syms.AddSymbol("AAG3", 109);
    dna_syms.AddSymbol("AAT1", 110);
    dna_syms.AddSymbol("AAT2", 111);
    dna_syms.AddSymbol("AAT3", 112);
    dna_syms.AddSymbol("ACA1", 113);
    dna_syms.AddSymbol("ACA2", 114);
    dna_syms.AddSymbol("ACA3", 115);
    dna_syms.AddSymbol("ACC1", 116);
    dna_syms.AddSymbol("ACC2", 117);
    dna_syms.AddSymbol("ACC3", 118);
    dna_syms.AddSymbol("ACG1", 119);
    dna_syms.AddSymbol("ACG2", 120);
    dna_syms.AddSymbol("ACG3", 121);
    dna_syms.AddSymbol("ACT1", 122);
    dna_syms.AddSymbol("ACT2", 123);
    dna_syms.AddSymbol("ACT3", 124);
    dna_syms.AddSymbol("AGA1", 125);
    dna_syms.AddSymbol("AGA2", 126);
    dna_syms.AddSymbol("AGA3", 127);
    dna_syms.AddSymbol("AGC1", 128);
    dna_syms.AddSymbol("AGC2", 129);
    dna_syms.AddSymbol("AGC3", 130);
    dna_syms.AddSymbol("AGG1", 131);
    dna_syms.AddSymbol("AGG2", 132);
    dna_syms.AddSymbol("AGG3", 133);
    dna_syms.AddSymbol("AGT1", 134);
    dna_syms.AddSymbol("AGT2", 135);
    dna_syms.AddSymbol("AGT3", 136);
    dna_syms.AddSymbol("ATA1", 137);
    dna_syms.AddSymbol("ATA2", 138);
    dna_syms.AddSymbol("ATA3", 139);
    dna_syms.AddSymbol("ATC1", 140);
    dna_syms.AddSymbol("ATC2", 141);
    dna_syms.AddSymbol("ATC3", 142);
    dna_syms.AddSymbol("ATG1", 143);
    dna_syms.AddSymbol("ATG2", 144);
    dna_syms.AddSymbol("ATG3", 145);
    dna_syms.AddSymbol("ATT1", 146);
    dna_syms.AddSymbol("ATT2", 147);
    dna_syms.AddSymbol("ATT3", 148);
    dna_syms.AddSymbol("CAA1", 149);
    dna_syms.AddSymbol("CAA2", 150);
    dna_syms.AddSymbol("CAA3", 151);
    dna_syms.AddSymbol("CAC1", 152);
    dna_syms.AddSymbol("CAC2", 153);
    dna_syms.AddSymbol("CAC3", 154);
    dna_syms.AddSymbol("CAG1", 155);
    dna_syms.AddSymbol("CAG2", 156);
    dna_syms.AddSymbol("CAG3", 157);
    dna_syms.AddSymbol("CAT1", 158);
    dna_syms.AddSymbol("CAT2", 159);
    dna_syms.AddSymbol("CAT3", 160);
    dna_syms.AddSymbol("CCA1", 161);
    dna_syms.AddSymbol("CCA2", 162);
    dna_syms.AddSymbol("CCA3", 163);
    dna_syms.AddSymbol("CCC1", 164);
    dna_syms.AddSymbol("CCC2", 165);
    dna_syms.AddSymbol("CCC3", 166);
    dna_syms.AddSymbol("CCG1", 167);
    dna_syms.AddSymbol("CCG2", 168);
    dna_syms.AddSymbol("CCG3", 169);
    dna_syms.AddSymbol("CCT1", 170);
    dna_syms.AddSymbol("CCT2", 171);
    dna_syms.AddSymbol("CCT3", 172);
    dna_syms.AddSymbol("CGA1", 173);
    dna_syms.AddSymbol("CGA2", 174);
    dna_syms.AddSymbol("CGA3", 175);
    dna_syms.AddSymbol("CGC1", 176);
    dna_syms.AddSymbol("CGC2", 177);
    dna_syms.AddSymbol("CGC3", 178);
    dna_syms.AddSymbol("CGG1", 179);
    dna_syms.AddSymbol("CGG2", 180);
    dna_syms.AddSymbol("CGG3", 181);
    dna_syms.AddSymbol("CGT1", 182);
    dna_syms.AddSymbol("CGT2", 183);
    dna_syms.AddSymbol("CGT3", 184);
    dna_syms.AddSymbol("CTA1", 185);
    dna_syms.AddSymbol("CTA2", 186);
    dna_syms.AddSymbol("CTA3", 187);
    dna_syms.AddSymbol("CTC1", 188);
    dna_syms.AddSymbol("CTC2", 189);
    dna_syms.AddSymbol("CTC3", 190);
    dna_syms.AddSymbol("CTG1", 191);
    dna_syms.AddSymbol("CTG2", 192);
    dna_syms.AddSymbol("CTG3", 193);
    dna_syms.AddSymbol("CTT1", 194);
    dna_syms.AddSymbol("CTT2", 195);
    dna_syms.AddSymbol("CTT3", 196);
    dna_syms.AddSymbol("GAA1", 197);
    dna_syms.AddSymbol("GAA2", 198);
    dna_syms.AddSymbol("GAA3", 199);
    dna_syms.AddSymbol("GAC1", 200);
    dna_syms.AddSymbol("GAC2", 201);
    dna_syms.AddSymbol("GAC3", 202);
    dna_syms.AddSymbol("GAG1", 203);
    dna_syms.AddSymbol("GAG2", 204);
    dna_syms.AddSymbol("GAG3", 205);
    dna_syms.AddSymbol("GAT1", 206);
    dna_syms.AddSymbol("GAT2", 207);
    dna_syms.AddSymbol("GAT3", 208);
    dna_syms.AddSymbol("GCA1", 209);
    dna_syms.AddSymbol("GCA2", 210);
    dna_syms.AddSymbol("GCA3", 211);
    dna_syms.AddSymbol("GCC1", 212);
    dna_syms.AddSymbol("GCC2", 213);
    dna_syms.AddSymbol("GCC3", 214);
    dna_syms.AddSymbol("GCG1", 215);
    dna_syms.AddSymbol("GCG2", 216);
    dna_syms.AddSymbol("GCG3", 217);
    dna_syms.AddSymbol("GCT1", 218);
    dna_syms.AddSymbol("GCT2", 219);
    dna_syms.AddSymbol("GCT3", 220);
    dna_syms.AddSymbol("GGA1", 221);
    dna_syms.AddSymbol("GGA2", 222);
    dna_syms.AddSymbol("GGA3", 223);
    dna_syms.AddSymbol("GGC1", 224);
    dna_syms.AddSymbol("GGC2", 225);
    dna_syms.AddSymbol("GGC3", 226);
    dna_syms.AddSymbol("GGG1", 227);
    dna_syms.AddSymbol("GGG2", 228);
    dna_syms.AddSymbol("GGG3", 229);
    dna_syms.AddSymbol("GGT1", 230);
    dna_syms.AddSymbol("GGT2", 231);
    dna_syms.AddSymbol("GGT3", 232);
    dna_syms.AddSymbol("GTA1", 233);
    dna_syms.AddSymbol("GTA2", 234);
    dna_syms.AddSymbol("GTA3", 235);
    dna_syms.AddSymbol("GTC1", 236);
    dna_syms.AddSymbol("GTC2", 237);
    dna_syms.AddSymbol("GTC3", 238);
    dna_syms.AddSymbol("GTG1", 239);
    dna_syms.AddSymbol("GTG2", 240);
    dna_syms.AddSymbol("GTG3", 241);
    dna_syms.AddSymbol("GTT1", 242);
    dna_syms.AddSymbol("GTT2", 243);
    dna_syms.AddSymbol("GTT3", 244);
    dna_syms.AddSymbol("TAA1", 245);
    dna_syms.AddSymbol("TAA2", 246);
    dna_syms.AddSymbol("TAA3", 247);
    dna_syms.AddSymbol("TAC1", 248);
    dna_syms.AddSymbol("TAC2", 249);
    dna_syms.AddSymbol("TAC3", 250);
    dna_syms.AddSymbol("TAG1", 251);
    dna_syms.AddSymbol("TAG2", 252);
    dna_syms.AddSymbol("TAG3", 253);
    dna_syms.AddSymbol("TAT1", 254);
    dna_syms.AddSymbol("TAT2", 255);
    dna_syms.AddSymbol("TAT3", 256);
    dna_syms.AddSymbol("TCA1", 257);
    dna_syms.AddSymbol("TCA2", 258);
    dna_syms.AddSymbol("TCA3", 259);
    dna_syms.AddSymbol("TCC1", 260);
    dna_syms.AddSymbol("TCC2", 261);
    dna_syms.AddSymbol("TCC3", 262);
    dna_syms.AddSymbol("TCG1", 263);
    dna_syms.AddSymbol("TCG2", 264);
    dna_syms.AddSymbol("TCG3", 265);
    dna_syms.AddSymbol("TCT1", 266);
    dna_syms.AddSymbol("TCT2", 267);
    dna_syms.AddSymbol("TCT3", 268);
    dna_syms.AddSymbol("TGA1", 269);
    dna_syms.AddSymbol("TGA2", 270);
    dna_syms.AddSymbol("TGA3", 271);
    dna_syms.AddSymbol("TGC1", 272);
    dna_syms.AddSymbol("TGC2", 273);
    dna_syms.AddSymbol("TGC3", 274);
    dna_syms.AddSymbol("TGG1", 275);
    dna_syms.AddSymbol("TGG2", 276);
    dna_syms.AddSymbol("TGG3", 277);
    dna_syms.AddSymbol("TGT1", 278);
    dna_syms.AddSymbol("TGT2", 279);
    dna_syms.AddSymbol("TGT3", 280);
    dna_syms.AddSymbol("TTA1", 281);
    dna_syms.AddSymbol("TTA2", 282);
    dna_syms.AddSymbol("TTA3", 283);
    dna_syms.AddSymbol("TTC1", 284);
    dna_syms.AddSymbol("TTC2", 285);
    dna_syms.AddSymbol("TTC3", 286);
    dna_syms.AddSymbol("TTG1", 287);
    dna_syms.AddSymbol("TTG2", 288);
    dna_syms.AddSymbol("TTG3", 289);
    dna_syms.AddSymbol("TTT1", 290);
    dna_syms.AddSymbol("TTT2", 291);
    dna_syms.AddSymbol("TTT3", 292);
}

#endif
