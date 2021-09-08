# coati

Branch | Status
------ | ------
CartwrightLab/master | [![github-actions](https://github.com/CartwrightLab/coati/actions/workflows/cmake.yml/badge.svg?branch=master)](https://github.com/CartwrightLab/coati/actions/workflows/cmake.yml) [![codecov](https://codecov.io/gh/CartwrightLab/coati/branch/master/graph/badge.svg)](https://codecov.io/gh/CartwrightLab/coati)
jgarciamesa/coati | [![github-actions](https://github.com/jgarciamesa/coati/actions/workflows/cmake.yml/badge.svg?branch=main)](https://github.com/jgarciamesa/coati/actions/workflows/meson.yml) [![codecov](https://codecov.io/gh/jgarciamesa/coati/branch/main/graph/badge.svg)](https://codecov.io/gh/jgarciamesa/coati)

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
meson setup builddir
meson compile -C builddir -Dbuildtype=release
```

### Global Install (requires root access)
```
tar -xvzf coati*.tar.gz
meson setup builddir
meson compile -C builddir -Dbuildtype=release
meson install -C builddir
```

## `alignpair`

Pairwise alignment of two nucleotide sequences. Example input files can be found in `fasta` directory.
```
Usage:	coati alignpair file.fasta [options]

Allowed options:
  -h [ --help ]                   Display this message
  -f [ --fasta ] arg              fasta file path
  -m [ --model ] arg (=m-coati)   substitution model: coati, m-coati (default),
                                  dna, ecm, m-ecm
  -w [ --weight ] arg             Write alignment score to file
  -o [ --output ] arg             Alignment output file
  -s [ --score ]                  Calculate alignment score using m-coati or
                                  m-ecm models
  -r [ --rate ] arg               Substitution rate matrix (CSV)
  -t [ --evo-time ] arg (=0.0133) Evolutionary time or branch length
```

### Sample runs:

```
#Align file examle-001.fasta with m-coati model (default) and output in PHY format (default)
coati alignpair fasta/example-001.fasta

# Align file example-002.fasta with m-ecm model and output in fasta format
coati alignpair fasta/example-002.fasta -m m-ecm -o example-002.fasta

# Align file example-003.fasta with ecm model, PHY output, and save alignment weight to w.out
coati alignpair fasta/example-003.fasta -m ecm -w w.out
```
