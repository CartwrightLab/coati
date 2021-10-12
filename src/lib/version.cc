/*
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
*/

#include <doctest/doctest.h>

#include <coati/coati.hpp>

// do some version number sanity checks
static_assert(COATI_VERSION_MAJOR >= 0 && COATI_VERSION_MAJOR < 1000,  // NOLINT
              "COATI major version must be less than 1000.");
static_assert(COATI_VERSION_MINOR >= 0 && COATI_VERSION_MINOR < 1000,  // NOLINT
              "COATI minor version must be less than 1000.");
static_assert(COATI_VERSION_PATCH >= 0 &&
                  COATI_VERSION_PATCH < 10000,  // NOLINT
              "COATI patch version must be less than 1000.");

/// @private
bool coati::version_number_check_equal(int version_int) {
    return version_int == COATI_VERSION_INTEGER;
}

/// @private
TEST_CASE("version_number_check_equal") {
    CHECK(coati::version_number_check_equal(COATI_VERSION_INTEGER) == true);
    CHECK(coati::version_number_check_equal(-1) == false);
}

/// @private
int coati::version_integer() { return COATI_VERSION_INTEGER; }

/// @private
TEST_CASE("version_integer") {
    CHECK(coati::version_integer() == COATI_VERSION_INTEGER);
}
