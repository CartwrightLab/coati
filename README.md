# coati

Branch | Status
------ | ------
CartwrightLab/master | [![github-actions](https://github.com/CartwrightLab/coati/actions/workflows/cmake.yml/badge.svg?branch=master)](https://github.com/CartwrightLab/coati/actions/workflows/cmake.yml) [![codecov](https://codecov.io/gh/CartwrightLab/coati/branch/master/graph/badge.svg)](https://codecov.io/gh/CartwrightLab/coati)
jgarciamesa/coati | [![github-actions](https://github.com/jgarciamesa/coati/actions/workflows/meson.yml/badge.svg?branch=main)](https://github.com/jgarciamesa/coati/actions/workflows/meson.yml) [![codecov](https://codecov.io/gh/jgarciamesa/coati/branch/main/graph/badge.svg)](https://codecov.io/gh/jgarciamesa/coati)

Codon-Aware Multiple Sequence Alignments

## Table of Contents
* [Installation](#installation)
* [alignpair](#alignpair)

## Installation

### Download
Source code for the most recent beta versions is available at <https://github.com/CartwrightLab/coati/archive/master.tar.gz>

### Dependencies

* Recent C++ compiler, supporting C++11 (e.g. gcc 4.8.1+ or clang 3.3+)
* Meson 0.59+
* Ninja 1.8.2+
* Boost 1.47+ <http://www.boost.org/>

### Compiling
```
tar -xvzf coati*.tar.gz
meson setup builddir --buildtype=release
meson compile -C builddir
```

### Global Install (requires root access)
```
tar -xvzf coati*.tar.gz
meson setup builddir --buildtype=release
meson compile -C builddir
meson install -C builddir
```

## `alignpair`

Pairwise alignment of two nucleotide sequences. Example input files can be found in `fasta` directory.
```
Usage: coati alignpair [OPTIONS] fasta

Positionals:
  fasta TEXT:FILE REQUIRED    Fasta file path

Options:
  -h,--help                   Print this help message and exit
  -m,--model TEXT             Substitution model
  -t,--time FLOAT:POSITIVE    Evolutionary time/branch length
  -l,--weight TEXT            Write alignment score to file
  -s,--score                  Score alignment
  -o,--output TEXT            Alignment output file
  -g,--gap-open FLOAT:POSITIVE
                              Gap opening score
  -e,--gap-extend FLOAT:POSITIVE
                              Gap extension score
  -w,--omega FLOAT:POSITIVE   Nonsynonymous-synonymous bias
  -p,--pi FLOAT x 4           Nucleotide frequencies (A C G T)
  -n,--gap-len UINT           Set gap unit size
```

### Sample runs:

```
#Align file examle-001.fasta with m-coati model (default) and output in PHY format (default)
coati alignpair fasta/example-001.fasta

# Align file example-002.fasta with m-ecm model and output in fasta format
coati alignpair fasta/example-002.fasta -m m-ecm -o example-002.fasta

# Align file example-003.fasta with ecm model, PHY output, and save alignment weight to w.out
coati alignpair fasta/example-003.fasta -m ecm -l w.out

# Align file example-003.fasta with m-coati model, specifying a branch length of 0.0133 and nucleotide frequencies A:0.3, C:0.2, C:0.2, T:0.3
coati alignpair fasta/example-003.fasta -t 0.0133 -p 0.3 0.2 0.2 0.3
```
