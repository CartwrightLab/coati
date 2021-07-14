# Copyright (c) 2020 Reed A. Cartwright <reed@cartwright.ht>
# Copyright (c) 2021 Juan J. Garcia Mesa <juanjosegarciamesa@gmail.com>
#
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

set(coati_sources version.cc mutation_coati.cc utils.cc align.cc tree.cc profile_aln.cc insertions.cc mutation_ecm.cc mutation_fst.cc gotoh.cc)
set(coati_headers coati.hpp mutation_coati.hpp utils.hpp align.hpp tree.hpp profile_aln.hpp dna_syms.hpp insertions.hpp mutation_ecm.hpp mutation_fst.hpp gotoh.hpp mg94p.tcc mg94q.tcc)

#####################################################################
# libcoati library
foreach(source IN LISTS coati_headers)
    list(APPEND coati_headers_ "${CMAKE_SOURCE_DIR}/src/include/coati/${source}")
endforeach()
set(coati_headers ${coati_headers_})

add_library(libcoati STATIC ${coati_sources} ${coati_headers})
add_dependencies(libcoati configure-version.h)
target_compile_features(libcoati PUBLIC cxx_std_17)
target_include_directories(libcoati PUBLIC "${CMAKE_SOURCE_DIR}/src/include")
target_include_directories(libcoati PUBLIC "${CMAKE_CURRENT_BINARY_DIR}/..")
target_compile_definitions(libcoati PRIVATE DOCTEST_CONFIG_DISABLE)

target_link_libraries(libcoati PRIVATE doctest::doctest)
target_link_libraries(libcoati PRIVATE minionrng::minionrng)
target_link_libraries(libcoati PRIVATE Boost::program_options)
target_link_libraries(libcoati PRIVATE Eigen3::Eigen)
target_link_libraries(libcoati PRIVATE FSTLIB::fst)

set_target_properties(libcoati PROPERTIES OUTPUT_NAME coati)

add_library(libcoati-doctest EXCLUDE_FROM_ALL OBJECT ${coati_sources})
target_compile_features(libcoati-doctest PUBLIC cxx_std_17)
target_link_libraries(libcoati-doctest PUBLIC
     $<TARGET_PROPERTY:libcoati,LINK_LIBRARIES>
)
target_include_directories(libcoati-doctest PRIVATE
    $<TARGET_PROPERTY:libcoati,INCLUDE_DIRECTORIES>
)
target_include_directories(libcoati-doctest INTERFACE
    $<TARGET_PROPERTY:libcoati,INTERFACE_INCLUDE_DIRECTORIES>
)

add_dependencies(libcoati-doctest configure-version.h)
target_link_libraries(libcoati-doctest PUBLIC devopt_coverage)

clang_format_target(libcoati)
clang_tidy_target(libcoati)
