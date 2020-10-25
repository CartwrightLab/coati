# Developer Documentation

## Branches

 - master is the development branch
 - stable is the releases branch
 - stable branch is tagged for releases

## Versions

 - the master branch will only have major and minor versions defined
 - when a new release/stable branch is made from master, the patch version will be 0, and the minor version will be incremented on master
 - version compiled directly from master will be known as prereleases

## Targets to run checks and tests
 - `make test` runs the unit and integration tests via `ctest`.
 - `make check_format` tests source code is formatted properly via [`clang-format`](https://clang.llvm.org/docs/ClangFormat.html). Use `make format` to fix format issues.
 - `make check_cppcheck` tests for bugs in the source code via [`Cppcheck`](http://cppcheck.sourceforge.net/).
    Use `make cppcheck` to print a bug report.
 - `make check_tidy` tests for bugs in the source code via [`clang-tidy`](https://clang.llvm.org/extra/clang-tidy/).
    Use `make tidy` to print a bug report.
 - `make check_all` runs all the above tests and is useful for testing code before submitting a pull request.
 
### Agave
Steps to set up COATi `only-gotoh` branch on the ASU Agave server:

```
module load cmake/latest
module load boost/1.61.0-gcc-9.2.0
module load eigen/3.3.7-gcc-stock

git clone https://github.com/jgarciamesa/coati.git --branch only-gotoh
cd coati/build
cmake -DFSTLIB_ROOT=~/tools/openfst-1.7.9/src/ -DEigen3_DIR=/packages/7x/eigen/3.3.7-gcc-stock/share/eigen3/cmake/ -DCMAKE_BUILD_TYPE=Release ..
make -j4
```

Example input files can be found in `fasta` directory. Sample runs:

```
#Align file ENSG00000004534.fasta with m-coati model and output in PHY format
./build/src/coati-alignpair fasta/ENSG00000004534.fasta

# Align file ENSG00000003249.fasta with m-ecm model and output in fasta format
./build/src/coati-alignpair fasta/ENSG00000003249.fasta -m m-ecm -o ENSG00000003249.fasta

# Align file ENSG00000005175.fasta with ecm model, PHY output, and save weight to w.out
./build/src/coati-alignpair fasta/ENSG00000005175.fasta -m ecm -w w.out

# Calculate the score of the second alignment with m-coati model
./build/src/coati-alignpair ENSG00000003249.fasta -s
```
