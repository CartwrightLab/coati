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
 