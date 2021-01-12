# Copyright (c) 2020 Reed A. Cartwright <reed@cartwright.ht>
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

include(CheckIncludeFile)
include(CheckFunctionExists)
include(CheckLibraryExists)

add_subdirectory(lib)

# Setup a target that will build version.h.in -> version.h
set(configure_script_code "
include(\"${CMAKE_CURRENT_BINARY_DIR}/../${CMAKE_PROJECT_NAME}.version\")\n
")
add_configured_file(version.h
    INPUT "${CMAKE_CURRENT_SOURCE_DIR}/version.h.in"
    OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/version.h"
    EXTRA_DEPS "${CMAKE_CURRENT_BINARY_DIR}/../${CMAKE_PROJECT_NAME}.version"
    VARIABLES ${CMAKE_PROJECT_NAME}_SEMVER_VERSION
              ${CMAKE_PROJECT_NAME}_SEMVER_INTEGER
              ${CMAKE_PROJECT_NAME}_SEMVER_MAJOR
              ${CMAKE_PROJECT_NAME}_SEMVER_MINOR
              ${CMAKE_PROJECT_NAME}_SEMVER_PATCH
)
add_dependencies(configure-version.h configure-${CMAKE_PROJECT_NAME}.version)
unset(configure_script_code)

# coati-version
add_executable(coati-version coati-version.cc)
target_link_libraries(coati-version libcoati)
install(TARGETS coati-version RUNTIME DESTINATION ${CMAKE_INSTALL_LIBEXECDIR})

clang_format_target(coati-version)
clang_tidy_target(coati-version)

# # coati-genseed
# add_executable(coati-genseed coati-genseed.cc)
# target_link_libraries(coati-genseed libcoati)
# target_link_libraries(coati-genseed minionrng::minionrng)
# install(TARGETS coati-genseed RUNTIME DESTINATION ${CMAKE_INSTALL_LIBEXECDIR})
#
# clang_format_target(coati-genseed)
# clang_tidy_target(coati-genseed)

# coati-alignpair
add_executable(coati-alignpair coati-alignpair.cc)
target_link_libraries(coati-alignpair libcoati)
target_link_libraries(coati-alignpair Boost::program_options)
target_link_libraries(coati-alignpair Boost::filesystem)
target_link_libraries(coati-alignpair Eigen3::Eigen)
target_link_libraries(coati-alignpair FSTLIB::fst)
install(TARGETS coati-alignpair RUNTIME DESTINATION ${CMAKE_INSTALL_LIBEXECDIR})

# coati main script
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/coati.cc.in"
  "${CMAKE_CURRENT_BINARY_DIR}/coati.cc"
  @ONLY)
add_executable(coati "${CMAKE_CURRENT_BINARY_DIR}/coati.cc")
target_link_libraries(coati libcoati)
install(TARGETS coati RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

###########################################################
# Unit Test support
add_executable(coati-doctest EXCLUDE_FROM_ALL coati-doctest.cc)

target_link_libraries(coati-doctest libcoati-doctest doctest::doctest)

doctest_discover_tests(coati-doctest PROPERTIES "FIXTURES_REQUIRED;doctest_fixture")

add_test(NAME [=["[pretest] build_coati-doctest"]=] COMMAND "${CMAKE_COMMAND}"
  --build "${CMAKE_BINARY_DIR}"
  --config "$<CONFIG>"
  --target coati-doctest
)
set_tests_properties([=["[pretest] build_coati-doctest"]=] PROPERTIES FIXTURES_SETUP doctest_fixture)

add_dependencies(pretest coati-doctest)
