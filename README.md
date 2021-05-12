# coati
Codon-Aware Multiple Sequence Alignments

## Table of Contents
* [Installation](#installation)
* [alignpair](#alignpair)

## Installation

### Download
Source code for the most recent beta versions is available at <https://github.com/CartwrightLab/coati/archive/master.tar.gz>

### Dependencies

* Recent C++ compiler, supporting C++11 (e.g. gcc 4.8.1+ or clang 3.3+)
* CMake 3.12+ when compiling <http://www.cmake.org/download/#latest>
* Boost 1.47+ <http://www.boost.org/>
* Eigen 3.3+ <http://eigen.tuxfamily.org/>
* OpenFST 1.8.0+ <http://openfst.org/twiki/bin/view/FST/FstDownload>


### Compiling
```
tar -xvzf coati*.tar.gz
cd coati*/build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j4
```

#### Compiling Flags
* `-DFSTLIB_ROOT` specify include directory for OpenFST library

### Global Install (requires root access)
```
cd coati*/build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j4
make install
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
