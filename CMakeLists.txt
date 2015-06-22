########################################################
#
#    Copyright (c) 2014-2015
#      SMASH Team
#
#    BSD 3-clause license
# 
######################################################### 

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

# Set a default value for CMAKE_BUILD_TYPE, otherwise we get something
# which is none of the options.
if(NOT CMAKE_BUILD_TYPE)
   # Change RelWithDebInfo below to Release for SMASH releases
   set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING
      "Choose the type of build, options are: Debug Release RelWithDebInfo Profiling."
      FORCE)
endif(NOT CMAKE_BUILD_TYPE)

# RelWithDebInfo is just like Release, it only adds -g to add debug symbols
# (useful for better debug and profiling information on Release builds)
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELEASE} -g" CACHE
   STRING "Flags used by the compiler during Release with Debug Info builds." FORCE)
set(CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO}" CACHE
   STRING "Flags used by the compiler during Release with Debug Info builds." FORCE)

# set up the profiling build type
set(CMAKE_CXX_FLAGS_PROFILING "-O3 -DNDEBUG -pg" CACHE STRING "Flags used by the compiler during profile builds." FORCE)
set(CMAKE_C_FLAGS_PROFILING "-O3 -DNDEBUG -pg" CACHE STRING "Flags used by the compiler during profile builds." FORCE)
set(CMAKE_EXE_LINKER_FLAGS_PROFILING "${CMAKE_EXE_LINKER_FLAGS_RELEASE} -pg" CACHE STRING "Flags used by the compiler during profile builds." FORCE)
mark_as_advanced(CMAKE_CXX_FLAGS_PROFILING CMAKE_C_FLAGS_PROFILING CMAKE_EXE_LINKER_FLAGS_RELEASE)

# Set the relevant generic compiler flags (optimisation + warnings)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W -Wall -Wextra -Wmissing-declarations -Wpointer-arith -Wshadow -Wuninitialized -Winit-self -Wfloat-equal -Wundef -Wcast-align -Wformat=2 -Wold-style-cast -std=c++11 -mfpmath=sse")

# Tell the compiler to ignore errno setting of math functions. This can help the compiler
# considerably in optimizing mathematical expressions.
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-math-errno")

if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
   # Clang-specific compiler flags
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-missing-braces")

   option(CLANG_USE_LIBC++ "If turned on clang will explicitly be asked to use libc++ (otherwise it uses the system default)" OFF)
   if(CLANG_USE_LIBC++)
      include(cmake/AddCompilerFlag.cmake)
      AddCompilerFlag(-stdlib=libc++ CXX_FLAGS CMAKE_CXX_FLAGS CXX_RESULT _use_libcxx)
      if(_use_libcxx AND "${CMAKE_SYSTEM_NAME}" STREQUAL "Linux")
         link_libraries(c++abi supc++)
      endif()
   endif()
endif()

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

# have binary in the build directory
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})

include_directories(SYSTEM "${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/einhard")

# enable standard CTest
# tests are enabled/disabled with the BUILD_TESTING switch to cmake
# uncomment the next line to default BUILD_TESTING to OFF
#option(BUILD_TESTING "Build the testing tree." OFF)
include(CTest)

# subdirectories where the code is
add_subdirectory(3rdparty)
add_subdirectory(src)
add_subdirectory(doc)
