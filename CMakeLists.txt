# Set CMake minimum required version
cmake_minimum_required(VERSION 3.1.0 FATAL_ERROR)

# Set default build type
SET(CMAKE_BUILD_TYPE_INIT "RelWithDebInfo")

SET(CMAKE_CXX_STANDARD 11)

# Set project name
project(coati)

# Find Boost installation with library program_options
find_package(Boost 1.47.0 REQUIRED COMPONENTS program_options)

#Find Eigen library
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

# Find FST library
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lfst" CACHE STRING "Set C++ Compiler Flags" FORCE)

# Check if Boost was found
if(Boost_FOUND)
	message("-- Boost found")
else()
	message (FATAL_ERROR "Cannot find Boost")
endif()

# Check if Eigen was found
if(TARGET Eigen3::Eigen)
	message("-- Eigen found")
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
	PRIVATE
		${PROJECT_SOURCE_DIR}/src/include/coati
)

# Link against Boost libraries
target_link_libraries(coati
	PRIVATE
		Boost::program_options
		Eigen3::Eigen
)
