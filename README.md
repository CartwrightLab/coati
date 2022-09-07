# coati

Branch | Status
------ | ------
CartwrightLab/master | [![github-actions](https://github.com/CartwrightLab/coati/actions/workflows/meson.yml/badge.svg?branch=master)](https://github.com/CartwrightLab/coati/actions/workflows/meson.yml) [![codecov](https://codecov.io/gh/CartwrightLab/coati/branch/master/graph/badge.svg)](https://codecov.io/gh/CartwrightLab/coati)
jgarciamesa/coati | [![github-actions](https://github.com/jgarciamesa/coati/actions/workflows/meson.yml/badge.svg?branch=main)](https://github.com/jgarciamesa/coati/actions/workflows/meson.yml) [![codecov](https://codecov.io/gh/jgarciamesa/coati/branch/main/graph/badge.svg)](https://codecov.io/gh/jgarciamesa/coati)

Codon-Aware Multiple Sequence Alignments

## Table of Contents
* [Installation](#installation)
* [alignpair](#alignpair)
* [format](#format)
* [msa](#msa)

## Installation

### Download
Source code for the most recent beta versions is available at <https://github.com/CartwrightLab/coati/archive/master.tar.gz>

### Dependencies

* Recent C++ compiler, supporting C++17 (e.g. gcc 8+ or clang 5+)
* Meson 0.59+ <https://mesonbuild.com/>
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

Syntax for output files is `[format:]filename.extension` where format is
optional (indicated by `[]`). Format must be one of "fa", "fasta" (FASTA format)
, "phy" (PHYLIP format), or "json". If no format is specified, then extension
must also be one of the above.

```
Usage: coati alignpair [OPTIONS] input

Positionals:
  input TEXT REQUIRED         Input file (FASTA/PHYLIP/JSON accepted)"

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
  -k,--gap-len UINT           Set gap unit size
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

## `format`

Convert between FASTA, PHYLIP, and JSON formats. PHYLIP format assums that all
sequences have equal length and JSON format is as follows:
```
{
    "data" : {
        "names" : ["Name_of_sequence1", "Name_of_sequence2"],
        "seqs" : ["AAAAAA", "AAAAAA"]
    }
}
```

Additionaly, `coati format` can adjust sequences aligned with `coati alignpair`
or `coati msa` to be used with downstream analyses and maintain our model
assumption that the reference (first sequence) is always in frame.
```
Usage: coati format [OPTIONS] input

Positionals:
  input TEXT REQUIRED         Input file (FASTA/PHYLIP/JSON accepted)

Options:
  -h,--help                   Print this help message and exit
  -o,--output TEXT            Alignment output file
  -p,--preserve-phase         Preserve phase
  -c,--padding TEXT Needs: --preserve-phase
                              Padding char to format preserve phase
```

## `msa`

Multiple sequence analyses (still in development).
```
Usage: coati msa [OPTIONS] input tree reference

Positionals:
  input TEXT REQUIRED         Input file (FASTA/PHYLIP/JSON accepted)
  tree TEXT:FILE REQUIRED     Newick phylogenetic tree
  reference TEXT REQUIRED     Reference sequence

Options:
  -h,--help                   Print this help message and exit
  -m,--model TEXT             Substitution model
  -o,--output TEXT            Alignment output file
  -g,--gap-open FLOAT:POSITIVE
                              Gap opening score
  -e,--gap-extend FLOAT:POSITIVE
                              Gap extension score
  -w,--omega FLOAT:POSITIVE   Nonsynonymous-synonymous bias
  -p,--pi FLOAT x 4           Nucleotide frequencies (A C G T)
  -k,--gap-len UINT           Set gap unit size
  ```
