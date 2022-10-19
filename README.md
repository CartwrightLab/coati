[![github-actions](https://github.com/CartwrightLab/coati/actions/workflows/meson.yml/badge.svg?branch=master)](https://github.com/CartwrightLab/coati/actions/workflows/meson.yml) [![codecov](https://codecov.io/gh/CartwrightLab/coati/branch/master/graph/badge.svg)](https://codecov.io/gh/CartwrightLab/coati)

# COATi - Codon-Aware Multiple Sequence Alignments

A statistical pairwise and multiple (in development) sequence aligner that is
robust to artifacts, incorporates codon models, and supports complex indels.

## Table of Contents
* [Installation](#installation)
* [`alignpair` - Pairwise alignment](#alignpair)
* [`sample` - Sample alignments](#sample)
* [`genseed` - Generate a random seed](#genseed)
* [`format` - Utility command to format sequences](#format)
* [`msa` - Multiple sequence aligner](#msa)
* [Input/output syntax details](#input-output-syntax)

## Installation

### Download
Source code for the most recent beta versions is available at <https://github.com/CartwrightLab/coati/archive/master.tar.gz>

### Dependencies

* Recent C++ compiler, supporting C++17 (e.g. gcc 8+ or clang 5+)
* [Meson 0.59+](https://mesonbuild.com/)
* [Ninja 1.8.2+](https://ninja-build.org/)
* [Boost 1.47+](http://www.boost.org/)

### Compiling
```bash
tar -xvzf coati*.tar.gz
meson setup builddir --buildtype=release
meson compile -C builddir
```

### Global Install (requires root access)
```bash
tar -xvzf coati*.tar.gz
meson setup builddir --buildtype=release
meson compile -C builddir
meson install -C builddir
```

## alignpair

Pairwise alignment of nucleotide sequences via composition of finite-state
transducers (FSTs) or using Gotoh's algorithm.
Substitution models available:

  - coati: Muse and Gaut codon model using finite-state transducers.
  - ecm: empirical codon model using FSTs.
  - dna: dna model resulting of marginalizing Muse and Gaut's model using FSTs.
  - m-coati: marginal Muse and Gaut using Gotoh's algorithm.
  - m-ecm: marginal empirical codon model using Gotoh's algorithm.

Indel model allows gaps of any length at any position.
Insertion always precede deletions when contiguous to eliminate equivalent
alignments.

```
coati alignpair - pairwise alignment of nucleotide sequences

Usage: coati alignpair [OPTIONS] input

Positionals:
  input TEXT REQUIRED              Input file (FASTA/PHYLIP/JSON accepted)

Options:
  -h,--help                        Print this help message and exit
  -l,--weight TEXT                 Write alignment score to file
  -s,--score                       Score input alignment and exit
  -o,--output TEXT                 Alignment output file


Model parameters:
  -m,--model TEXT Excludes: --sub  Substitution model (coati ecm dna m-coati m-ecm)
  -t,--time FLOAT:POSITIVE         Evolutionary time/branch length
  -g,--gap-open FLOAT:POSITIVE     Gap opening score
  -e,--gap-extend FLOAT:POSITIVE   Gap extension score
  -k,--gap-len UINT                Gap unit length


Advanced options:
  --sub TEXT Excludes: --model     File with branch lengths and codon subst matrix
  -r,--ref TEXT Excludes: --rev-ref
                                   Name of reference sequence (default: 1st seq)
  -v,--rev-ref Excludes: --ref     Use 2nd seq as reference (default: 1st seq)
  -w,--omega FLOAT:POSITIVE        Nonsynonymous-synonymous bias
  -p,--pi FLOAT x 4                Nucleotide frequencies (A C G T)
  -x,--sigma FLOAT x 6             GTR sigma parameters (AC AG AT CG CT GT)

```

### Sample runs

```bash
#Align file example-001.fasta with m-coati model (default) and output in JSON
format (default)
coati alignpair sampledata/example-001.fasta

# Align file example-002.fasta with m-ecm model and output in fasta format
coati alignpair sampledata/example-002.fasta -m m-ecm -o example-002.fasta

# Align file example-003.fasta with ecm model, PHY output, and save alignment weight to w.out
coati alignpair sampledata/example-003.fasta -m ecm -l w.out

# Align file example-003.fasta with m-coati model, specifying a branch length of 0.0133 and underlying nucleotide frequencies A:0.3, C:0.2, C:0.2, T:0.3
coati alignpair sampledata/example-003.fasta -t 0.0133 -p 0.3 0.2 0.2 0.3
```

## sample

Align a pair of sequences and get alignments from Gotoh's dynamic programming
matrices with sampling (not necessarily best alignment).
Output is always in json format ([see input output syntax](#input-output-syntax))
and provides a weight and log weight of each pairwise alignment sampled.

```
coati sample - align two sequences and sample alignments

Usage: coati sample [OPTIONS] input

Positionals:
  input TEXT REQUIRED              Input file (FASTA/PHYLIP/JSON accepted)

Options:
  -h,--help                        Print this help message and exit
  -o,--output TEXT                 Alignment output file
  -n,--sample-size UINT            Sample size
  -s,--seed TEXT ...               Space separated list of seed(s) used for sampling


Model parameters:
  -t,--time FLOAT:POSITIVE         Evolutionary time/branch length
  -m,--model TEXT Excludes: --sub  Substitution model (coati ecm dna m-coati m-ecm)
  -g,--gap-open FLOAT:POSITIVE     Gap opening score
  -e,--gap-extend FLOAT:POSITIVE   Gap extension score
  -k,--gap-len UINT                Gap unit length


Advanced options:
  --sub TEXT Excludes: --model     File with branch lengths and codon subst matrix
  -w,--omega FLOAT:POSITIVE        Nonsynonymous-synonymous bias
  -p,--pi FLOAT x 4                Nucleotide frequencies (A C G T)
  -x,--sigma FLOAT x 6             GTR sigma parameters (AC AG AT CG CT GT)
```

### Sample runs

```bash
# Sample 100 alignments from example-003.fasta and save to `ex3_100.json`.
coati sample sampledata/example-003.fasta -n 100 -o ex3_100.json

# Sample 50 alignments from example-003.fasta specifying seed
coati sample sampledata/example-003.fasta -n 50 -s random42
```

## genseed

Generate a random seed. Usage:

```bash
coati genseed
```

## format

Convert between formats, extract and/or reorder sequences.
Accepted format are FASTA, PHYLIP, and JSON formats.
Output can be piped to other commands.

Additionaly, `coati format` can adjust sequences aligned with `coati alignpair`
or `coati msa` to be used with downstream analyses and maintain our model
assumption that the reference is always in frame.
```
coati format - convert between formats, extract and/or reoder sequences

Usage: coati format [OPTIONS] input

Positionals:
  input TEXT REQUIRED              Input file (FASTA/PHYLIP/JSON accepted)

Options:
  -h,--help                        Print this help message and exit
  -o,--output TEXT                 Alignment output file
  -p,--preserve-phase              Preserve phase
  -c,--padding TEXT Needs: --preserve-phase
                                   Padding char to format preserve phase
  -s,--cut-seqs TEXT ... Excludes: --cut-pos
                                   Name of sequences to extract
  -x,--cut-pos UINT ... Excludes: --cut-seqs
                                   Position of sequences to extract (1 based)

```

### Sample runs

```bash
# Extract first two sequences and align them with alignpair
coati format sampledata/example-msa-001.fasta -x 1 2 | coati alignpair json:-

# Extract sequences C B D (in that order) and convert to PHYLIP format
coati format sampledata/example-msa-002.fasta -s C B D -o ex-msa-2.phy

# Align then preserve phase using '?' character
coati alignpair sampledata/example-003.fasta | coati format json:- -p -c '?'
```

## msa

Multiple sequence alignment (still in development).

```
coati msa - multiple sequence alignment of nucleotide sequences

Usage: coati msa [OPTIONS] input tree reference

Positionals:
  input TEXT REQUIRED              Input file (FASTA/PHYLIP/JSON accepted)
  tree TEXT:FILE REQUIRED          Newick phylogenetic tree
  reference TEXT REQUIRED          Name of reference sequence

Options:
  -h,--help                        Print this help message and exit
  -o,--output TEXT                 Alignment output file


Model parameters:
  -m,--model TEXT                  Substitution model (coati ecm dna m-coati m-ecm)
  -g,--gap-open FLOAT:POSITIVE     Gap opening score
  -e,--gap-extend FLOAT:POSITIVE   Gap extension score
  -k,--gap-len UINT                Gap unit length


Advanced options:
  -w,--omega FLOAT:POSITIVE        Nonsynonymous-synonymous bias
  -p,--pi FLOAT x 4                Nucleotide frequencies (A C G T)
  -x,--sigma FLOAT x 6             GTR sigma parameters (AC AG AT CG CT GT)

  ```

## Input output syntax

Syntax for input and output files is `[format:]filename.extension` where format
is optional (indicated by `[]`).
Format must be one of "fa", "fasta" (FASTA format), "phy" (PHYLIP format), or
"json" (JSON format).
If no format is specified, then extension must also be one of the above.
Alternatively, input can be of the form `format:-` or `format:` and COATi will
read from the command line (std::cin);
This allows the different COATi commands to work in a pipeline.

### JSON output

The JSON input and output (except for `sample`) JSON format is as follows:
```json
{
  "data": {
    "names": ["A","B","C","D","E"],
    "seqs":["TCA--TCG","TCA-GTCG","T-A--TCG","TCAC-TCG","TCA--TC-"]
  }
}
```

The output from `sample` is always in the following JSON format:

```json
[
  {
    "aln": {
      "1": "CTCTGGATAGTG",
      "2": "CT----ATAGTG"
    },
    "weight": 0.994785,
    "log_weight": -0.00522821
  },
  {
    "aln": {
      "1": "CTCTGGATAGTG",
      "2": "CT----ATAGTG"
    },
    "weight": 0.994785,
    "log_weight": -0.00522821
  },
  {
    "aln": {
      "1": "CTCTGGATAGTG",
      "2": "CT----ATAGTG"
    },
    "weight": 0.994785,
    "log_weight": -0.00522821
  }
]

```

### PHYLIP

PHYLIP format limits names to 10 characters and uses interleaved format:

```
3 110
A         CCTACTCACTAGGAACAATTCTAAGGTTAT--------------------
B         CCTACTCACGACGAACAAGTTTAAGGTTATATGACACTTCTCCGATTGTC
VeryLongNaCCTACTCACTAGGAACAATTCTAAGGTTAT--------------------

-----GTTCCACTTCACCAAGTGTCGCATGTCGTCCACAGGTGTTCATTG----GCTCAA
GCAAGG-------------------------CGTCTACAGGTGTTCATTGTTACG----A
-----GTTCCACTTCACCAAGTGTCGCATGTCGTCCACAGGTGTTCATTGTTACG----A
```
