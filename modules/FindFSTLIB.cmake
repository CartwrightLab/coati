# Copyright (c) 2020 Juan J. Garcia Mesa <juanjosegarciamesa@gmail.com>
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

# Inspired by https://github.com/Kitware/CMake/blob/master/Modules/FindZLIB.cmake

set(_FSTLIB_SEARCHES)

# Search FSTLIB_ROOT first if it is set
if(FSTLIB_ROOT)
  set(_FSTLIB_SEARCH_ROOT PATHS ${FSTLIB_ROOT} NO_DEFAULT_PATH)
  list(APPEND _FSTLIB_SEARCHES _FSTLIB_SEARCH_ROOT)
endif()

# Normal search
set(_FSTLIB_SEARCH_NORMAL
	PATHS "/usr/local/")
list(APPEND _FSTLIB_SEARCHES _FSTLIB_SEARCH_NORMAL)

set(FSTLIB_NAMES fst fstlib)

# Try each search configuration
foreach(search ${_FSTLIB_SEARCHES})
	find_path(FSTLIB_INCLUDE_DIR NAMES fst/fstlib.h ${${search}} PATH_SUFFIXES "include")
endforeach()

foreach(search ${_FSTLIB_SEARCHES})
	find_path(FSTLIB_VERSION_FILE NAMES libfst.la ${${search}} PATH_SUFFIXES "lib")
endforeach()



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


if(FSTLIB_INCLUDE_DIR AND EXISTS "${FSTLIB_VERSION_FILE}/libfst.la")
	file(STRINGS "${FSTLIB_VERSION_FILE}/libfst.la" FSTLIB_VERSION_MAJOR REGEX "^current=([0-9]+)")
	file(STRINGS "${FSTLIB_VERSION_FILE}/libfst.la" FSTLIB_VERSION_MINOR REGEX "^age=([0-9]+)")
	file(STRINGS "${FSTLIB_VERSION_FILE}/libfst.la" FSTLIB_VERSION_PATCH REGEX "^revision=([0-9]+)")

	string(REGEX REPLACE "current=" "" FSTLIB_VERSION_MAJOR "${FSTLIB_VERSION_MAJOR}")
	string(REGEX REPLACE "age=" "" FSTLIB_VERSION_MINOR  "${FSTLIB_VERSION_MINOR}")
	string(REGEX REPLACE "revision=" "" FSTLIB_VERSION_PATCH "${FSTLIB_VERSION_PATCH}")

	set(FSTLIB_VERSION "${FSTLIB_VERSION_MAJOR}.${FSTLIB_VERSION_MINOR}.${FSTLIB_VERSION_PATCH}")

	# check if version found is less than minumum version required
	if(${FSTLIB_VERSION} VERSION_LESS ${FSTLIB_FIND_VERSION})
		set(FSTLIB_VERSION_OK FALSE)
	else()
		set(FSTLIB_VERSION_OK TRUE)
	endif()

else()
	message(SEND_ERROR "FSTLIB not found.")
	return()
endif()

if(NOT FSTLIB_VERSION_OK)
	message(STATUS "FSTLIB version ${FSTLIB_VERSION} found in ${FSTLIB_INCLUDE_DIR}, "
               "but at least version ${FSTLIB_FIND_VERSION} is required")
else(NOT FSTLIB_VERSION_OK)
	set(FSTLIB_LIBRARY ${FSTLIB_INCLUDE_DIR}/fst/fstlib.h)
	message(STATUS "FSTLIB version ${FSTLIB_VERSION} found in ${FSTLIB_INCLUDE_DIR}")
endif(NOT FSTLIB_VERSION_OK)


include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FSTLIB FSTLIB_LIBRARY FSTLIB_INCLUDE_DIR)

if(FSTLIB_FOUND)
	set(FSTLIB_INCLUDE_DIRS ${FSTLIB_INCLUDE_DIR})

	if(NOT FSTLIB_LIBRARIES)
		set(FSTLIB_LIBRARIES ${FSTLIB_LIBRARY})
	endif()

	if(NOT TARGET FSTLIB::fst)
		add_library(FSTLIB::fst UNKNOWN IMPORTED)
		set_target_properties(FSTLIB::fst PROPERTIES INTERFACE_INCLUDE_DIRECTORIES ${FSTLIB_INCLUDE_DIR})
        set_property(TARGET FSTLIB::fst APPEND PROPERTY IMPORTED_LOCATION ${FSTLIB_LIBRARY})
    endif()
endif()
