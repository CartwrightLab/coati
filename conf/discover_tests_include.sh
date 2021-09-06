#!/bin/#

# find  includes for all tests
# $1 = directory

grep -R --no-filename '#include' "${1}" | sort | uniq | grep -v 'doctest'
