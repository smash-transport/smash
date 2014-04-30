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
set(VERSION_MAJOR 20140430)

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
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W -Wall -Wextra -Wmissing-declarations -Wpointer-arith -Wshadow -Wuninitialized -Winit-self -Wfloat-equal -Wundef -Wcast-align -Wformat=2 -std=c++11")

# check that C++11 code compiles
include(CheckCXXSourceCompiles)
CHECK_CXX_SOURCE_COMPILES("
#include <initializer_list>
struct A {};
struct B : public A { using A::A; };
template<typename T = void> int f() { return 1; }
int main() {
  auto l = []() { return f(); };
  return l() - 1;
}
" CXX11_COMPILES)
if(NOT CXX11_COMPILES)
   message(FATAL_ERROR "Your compiler did not manage to compile a simple C++11 program. Please get a newer C++ compiler.")
endif()

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
