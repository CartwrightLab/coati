# MIT License

# Copyright (c) 2020 Reed A. Cartwright <reed@cartwrig.ht>
# Copyright (c) 2021 Juan J. Garcia Mesa <juanjosegarciamesa@gmail.com>

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

add_library(doctest INTERFACE)
add_library(doctest::doctest ALIAS doctest)
target_include_directories(doctest INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/doctest")

add_library(minionrng STATIC minionrng/minion.cc)
add_library(minionrng::minionrng ALIAS minionrng)
target_include_directories(minionrng PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/minionrng")

set(fstlib_sources compat.cc encode.cc flags.cc fst.cc fst-types.cc mapped-file.cc properties.cc symbol-table.cc symbol-table-ops.cc util.cc weight.cc)
set(fstlib_headers "${CMAKE_CURRENT_SOURCE_DIR}/fstlib/fst/fstlib.h")

foreach(source IN LISTS fstlib_sources)
    list(APPEND fstlib_sources_ "${CMAKE_CURRENT_SOURCE_DIR}/fstlib/${source}")
endforeach()
set(fstlib_sources ${fstlib_sources_})

add_library(fst STATIC ${fstlib_sources} ${fstlib_headers})
add_library(FSTLIB::fst ALIAS fst)
target_include_directories(fst PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/fstlib")
