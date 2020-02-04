# Set CMake minimum required version
cmake_minimum_required(VERSION 3.1.0 FATAL_ERROR)

# Set default build type
SET(CMAKE_BUILD_TYPE_INIT "RelWithDebInfo")

# Set project name, version, description, and language
project(coati VERSION 0.0
	DESCRIPTION "COdon-aware Alignment Transducer - COATi"
	# HOMEPAGE_URL "https://github.com/"
	LANGUAGES CXX
)


CONFIGURE_FILE(coati_config.h.in coati_config.h)
SET(CMAKE_CXX_STANDARD 11)
SET(CMAKE_CXX_STANDARD_REQUIRED True)

# Find Boost installation with library program_options
find_package(Boost 1.47.0 REQUIRED COMPONENTS program_options unit_test_framework filesystem)

#Find Eigen library
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

# Find FST library
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lfst" CACHE STRING "Set C++ Compiler Flags" FORCE)

# Check if Eigen was found
if(TARGET Eigen3::Eigen)
	message(STATUS "Eigen3 v${EIGEN3_VERSION_STRING} found in ${EIGEN3_INCLUDE_DIR}")
else()
	message(FATAL_ERROR "Cannot find Eigen")
endif(TARGET Eigen3::Eigen)

# Enable unit and executable testing
ENABLE_TESTING()

# Create sources variable with a link to all cc files to compile
set(SOURCES
	src/main.cpp
	src/lib/mut_models.cpp
	src/lib/utils.cpp
)

# Add executable
add_executable(coati ${SOURCES})

# Set directories that should be included in the build command for this target
target_include_directories(coati
	PUBLIC
		"${PROJECT_BINARY_DIR}"
		"${PROJECT_SOURCE_DIR}/src/include/coati"
)

# Link against Boost libraries
target_link_libraries(coati
	PUBLIC
		Boost::program_options
		Boost::filesystem
		Eigen3::Eigen
)
