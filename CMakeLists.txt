# The name of our project
project(smash)

# Fail if cmake is called in the source directory
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
   message(FATAL_ERROR "You don't want to configure in the source directory!")
endif()

# Minimum cmake version this is tested on
cmake_minimum_required(VERSION 2.8)

# Tell cmake where to find our modules
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# The current version number
set(VERSION_MAJOR 20140124)

# Set the relevant generic g++ flags (optimisation + warnings)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W -Wall -Wextra -Wmissing-declarations -Wpointer-arith -Wshadow -Wuninitialized -Winit-self -Winline -Wfloat-equal -Wundef -Wcast-align -Wformat=2 -std=c++11")

# use C++ 11 standard
#set(CMAKE_CXX_FLAGS “${CMAKE_CXX_FLAGS} “)
#add_definitions()

# enable debug preprocessor variable
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG")

# have binary in the build directory
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})

# enable standard CTest
# tests are enabled/disabled with the BUILD_TESTING switch to cmake
# uncomment the next line to default BUILD_TESTING to OFF
#option(BUILD_TESTING "Build the testing tree." OFF)
include(CTest)

# tests are run with make test/ctest. Leave the check target in for a little while to tell where it went.
add_custom_target(check
   COMMENT "'make check' is deprecated. Use 'make test' or call ctest. The latter also supports running tests in parallel with -j"
   )

# subdirectories where the code is
add_subdirectory(src)
