cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(SMASH_as_a_library_example)

# Tell cmake where to find modules
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

add_executable(example example.cc)
add_executable(example_rate_equations example_rate_equations.cc)

# Set the relevant generic compiler flags (optimisation + warnings)
# Try to be close to what SMASH uses
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W -Wall -Wextra -Wmissing-declarations -Wpointer-arith -Wshadow -Wuninitialized -Winit-self -Wundef -Wcast-align -Wformat=2 -Wold-style-cast -Werror=switch")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -fno-math-errno -std=c++11")
message(STATUS "CXX_FLAGS = " ${CMAKE_CXX_FLAGS})


# TODO: This is a crutch, which sets big endian or little endian in the same
#       way that SMASH does. It is required to avoid sneaky errors with pdg codes
#       byte order. Ideally such thing should not be necessary -- smash has to set
#       everything up correctly on its side.
include(TestBigEndian)
TEST_BIG_ENDIAN(IS_BIG_ENDIAN)
if(IS_BIG_ENDIAN)
    message(STATUS "Big endian architecture detected.")
    add_definitions("-DBIG_ENDIAN_ARCHITECTURE")
else()
    message(STATUS "Little endian architecture detected.")
    add_definitions("-DLITTLE_ENDIAN_ARCHITECTURE")
endif()


find_package(SMASH)
if(${SMASH_FOUND})
  message("SMASH libraries = ${SMASH_LIBRARIES}")
  include_directories(${SMASH_INCLUDE_DIR})
  # SYSTEM supresses compiling warnings from external libraries
  target_link_libraries(example SYSTEM ${SMASH_LIBRARIES})
  target_link_libraries(example_rate_equations SYSTEM ${SMASH_LIBRARIES})
endif(${SMASH_FOUND})
