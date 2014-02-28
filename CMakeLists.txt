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
set(VERSION_MAJOR 20140228)

# Set the default for CMAKE_BUILD_TYPE to Release, otherwise we get something
# which is none of the options.
if(NOT CMAKE_BUILD_TYPE)
   # Change RelWithDebInfo below to Release for SMASH releases
   set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING
      "Choose the type of build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
endif(NOT CMAKE_BUILD_TYPE)

# RelWithDebInfo normally uses "-O2 -g -DNDEBUG" but we want -O3 and keep
# assertions enabled.
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g" CACHE
   STRING "Flags used by the compiler during Release with Debug Info builds." FORCE)
set(CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO}" CACHE
   STRING "Flags used by the compiler during Release with Debug Info builds." FORCE)

# Set the relevant generic g++ flags (optimisation + warnings)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
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

# subdirectories where the code is
add_subdirectory(src)
add_subdirectory(doc)
