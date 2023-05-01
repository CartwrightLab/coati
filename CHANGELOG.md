# COATi Release Notes and Changelog

## v1.0 (May 2023)

### Features
 - COATi has 6 available commands: `help`, `alignpair` - pairwise alignment,
   `msa` - multiple sequences alignment (in development), `sample` - align two
   sequences and sample alignments, `genseed` - generate a random seed, and
   `format` - convert format, extract and/or reorder sequences.
 - `coati alignpair` can now use five different models: tri-mg (triplet Muse and Gaut),
   mar-mg (marginal Muse and Gaut), tri-ecm (Empirical Codon Model), mar-ecm (marginal ECM),
   and dna.
 - Input and output formats now include FASTA, PHYLIP, and JSON.
 - COATi reads from stdin and writes to stdout (in JSON format) by default.
 - Users can now specify parameter values. For an up-to-date list of parameters
   run `coati-command --help`.
 - Ambiguous nucleotides on descendant sequence are supported.
 - Add doxygen documentation.

## Notes for Developers
 - FSTlib, CLI, and doctest dependencies are now included and built with COATi
   as static libraries.
 - Eigen, JSON library, and Google Benchmark are handled by Meson (via meson
   wrapDB).

## v0.1.0 (2020/02/18)

### Features
 - Initial Release
 - `coati version` returns the version information
 - `coati alignpair -f fasta.fas -m toycoati` aligns two sequences given in a FASTA formatted file using *toycoati* model.

### Notes for Developers
 - This release focused on setting up the project infrastructure including
   unit testing, continuous integration, code formatting, and code analysis.
