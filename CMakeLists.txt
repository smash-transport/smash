########################################################
#
#    Copyright (c) 2014-2015
#      SMASH Team
#
#    BSD 3-clause license
#
#########################################################

# The name of our project
project(smash CXX)

# Fail if cmake is called in the source directory
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
   message(FATAL_ERROR "You don't want to configure in the source directory!")
endif()

# Minimum cmake version this is tested on
cmake_minimum_required(VERSION 3.1)

if(NOT "${CMAKE_CXX_STANDARD}")
  set(CMAKE_CXX_STANDARD 11)
endif()
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Tell cmake where to find our modules
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# needed for clang-tidy
set(CMAKE_EXPORT_COMPILE_COMMANDS "ON")

# Set a default value for CMAKE_BUILD_TYPE, otherwise we get something
# which is none of the options.
if(NOT CMAKE_BUILD_TYPE)
   # Change RelWithDebInfo below to Release for SMASH releases
   set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING
      "Choose the type of build, options are: Debug Release RelWithDebInfo Profiling."
      FORCE)
endif(NOT CMAKE_BUILD_TYPE)


# add 3rd-party libraries (before setting compiler flags etc)
include_directories(SYSTEM "${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/Cuba-4.2")
include_directories(SYSTEM "${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/einhard")

find_package(yaml-cpp 0.6 QUIET NO_MODULE)
if (yaml-cpp_FOUND)
  # target doesn't set include directory for some reason
  get_filename_component(yaml_cpp_include_dir ${YAML_CPP_INCLUDE_DIR} ABSOLUTE)
  include_directories(SYSTEM ${yaml_cpp_include_dir})
else()
  include_directories(SYSTEM "${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/yaml-cpp-0.6.2/include")
endif()

add_subdirectory(3rdparty)

# Set compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mfpmath=sse")

# Disable DEBUG and TRACE output in Release builds
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -DNDEBUG" CACHE
   STRING "Flags used by the compiler during Release builds." FORCE)
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

# Set the relevant generic compiler flags (optimisations)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -fPIC")

# Tell the compiler to ignore errno setting of math functions. This can help the compiler
# considerably in optimizing mathematical expressions.
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-math-errno")

if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
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
# in order for this to work the policy on respecting the CXX_STANDARD for
# try_compile has to be changed first
# see https://cmake.org/cmake/help/v3.8/policy/CMP0067.html
cmake_policy(SET CMP0067 NEW)
set(CMAKE_REQUIRED_FLAGS "-std=c++11")
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

# enable standard CTest
# tests are enabled/disabled with the BUILD_TESTING switch to cmake
# uncomment the next line to default BUILD_TESTING to OFF
#option(BUILD_TESTING "Build the testing tree." OFF)
include(CTest)

# subdirectories where the code is
add_subdirectory(src)
add_subdirectory(doc)

install(DIRECTORY input/ DESTINATION share/smash)
