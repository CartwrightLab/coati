# Copyright (c) 2020 Juan J. Garcia Mesa <juanjosegarciamesa@gmail.com>
# Copyright (c) 2021 Reed A. Cartwright <cartwright@asu.edu>
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

# Try to find FstLib
# Once done this will define
# FSTLIB_FOUND - System has FstLib
# FSTLIB_INCLUDE_DIRS - FstLib include directories
# FSTLIB_LIBRARIES - The libraries needed to use FstLib

include(LibFindMacros)

find_path(FSTLIB_INCLUDE_DIR
	NAMES fst/fstlib.h
)

find_library(FST_LIBRARY
	NAMES fst
)

set(FSTLIB_PROCESS_INCLUDES FSTLIB_INCLUDE_DIR)
set(FSTLIB_PROCESS_LIBS FST_LIBRARY)
libfind_process(FSTLIB)

if(FSTLIB_FOUND)
	if(NOT TARGET FSTLIB::fst)
		add_library(FSTLIB::fst UNKNOWN IMPORTED)
		set_target_properties(FSTLIB::fst PROPERTIES INTERFACE_INCLUDE_DIRECTORIES ${FSTLIB_INCLUDE_DIR})
        set_property(TARGET FSTLIB::fst APPEND PROPERTY IMPORTED_LOCATION ${FSTLIB_LIBRARY})
    endif()
endif()
