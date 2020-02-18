# Copyright (c) 2020 Juan J. Garcia Mesa <jgarc111@asu.edu>
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

# Inspired by https://github.com/OPM/eigen3

# Find fstlib.h
find_path(FSTLIB_INCLUDE_DIR
	NAMES fstlib.h
	PATHS "/usr/local/include"
	PATH_SUFFIXES "fst"
	DOC "Include directory for Fstlibrary"
)

# Find libfst.la (contains version info)
find_path(FSTLIB_LIB_DIR
	NAMES libfst.la
	PATHS "/usr/local/lib"
	DOC "Library director for Fstlibrary"
)

# set minimum version to 17.0.0
if(NOT FSTLIB_FIND_VERSION)
	if(NOT FSTLIB_FIND_VERSION_MAJOR)
		set(FSTLIB_FIND_VERSION_MAJOR 17)
	endif(NOT FSTLIB_FIND_VERSION_MAJOR)
	if(NOT FSTLIB_FIND_VERSION_MINOR)
		set(FSTLIB_FIND_VERSION_MINOR 0)
	endif(NOT FSTLIB_FIND_VERSION_MINOR)
	if(NOT FSTLIB_FIND_VERSION_PATCH)
		set(FSTLIB_FIND_VERSION_PATCH 0)
	endif(NOT FSTLIB_FIND_VERSION_PATCH)

	set(FSTLIB_FIND_VERSION "${FSTLIB_FIND_VERSION_MAJOR}.${FSTLIB_FIND_VERSION_MINOR}.${FSTLIB_FIND_VERSION_PATCH}")

endif(NOT FSTLIB_FIND_VERSION)

if(FSTLIB_INCLUDE_DIR AND EXISTS "${FSTLIB_LIB_DIR}/libfst.la")
	file(STRINGS "${FSTLIB_LIB_DIR}/libfst.la" FSTLIB_VERSION_MAJOR REGEX "^current=([0-9]+)")
	file(STRINGS "${FSTLIB_LIB_DIR}/libfst.la" FSTLIB_VERSION_MINOR REGEX "^age=([0-9]+)")
	file(STRINGS "${FSTLIB_LIB_DIR}/libfst.la" FSTLIB_VERSION_PATCH REGEX "^revision=([0-9]+)")

	string(REGEX REPLACE "current=" "" FSTLIB_VERSION_MAJOR "${FSTLIB_VERSION_MAJOR}")
	string(REGEX REPLACE "age=" "" FSTLIB_VERSION_MINOR  "${FSTLIB_VERSION_MINOR}")
	string(REGEX REPLACE "revision=" "" FSTLIB_VERSION_PATCH "${FSTLIB_VERSION_PATCH}")

	set(FSTLIB_VERSION "${FSTLIB_VERSION_MAJOR}.${FSTLIB_VERSION_MINOR}.${FSTLIB_VERSION_PATCH}")

	# check if version found is less than minumum version required
	if(${FSTLIB_VERSION} VERSION_LESS ${FSTLIB_FIND_VERSION})
		set(FSTLIB_VERSION_OK FALSE)
	else(${FSTLIB_VERSION} VERSION_LESS ${FSTLIB_FIND_VERSION})
		set(FSTLIB_VERSION_OK TRUE)
	endif(${FSTLIB_VERSION} VERSION_LESS ${FSTLIB_FIND_VERSION})

	if(NOT FSTLIB_VERSION_OK)
		message(STATUS "Eigen3 version ${EIGEN3_VERSION} found in ${EIGEN3_INCLUDE_DIR}, "
	               "but at least version ${Eigen3_FIND_VERSION} is required")
	endif(NOT FSTLIB_VERSION_OK)
endif()


include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FSTLIB FSTLIB_INCLUDE_DIR)
