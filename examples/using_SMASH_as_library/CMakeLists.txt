cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(SMASH_as_a_library_example)

# Tell cmake where to find modules
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

add_executable(example example.cc)
add_executable(example_rate_equations example_rate_equations.cc)

# Set the relevant generic compiler flags (optimisation + warnings)
# Try to be close to what SMASH uses
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W -Wall -Wextra -Wmissing-declarations -Wpointer-arith -Wshadow -Wuninitialized -Winit-self -Wundef -Wcast-align -Wformat=2 -Wold-style-cast -Werror=switch")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-math-errno -std=c++11")
message(STATUS "CXX_FLAGS = " ${CMAKE_CXX_FLAGS})

find_package(SMASH)
if(${SMASH_FOUND})
  message("SMASH libraries = ${SMASH_LIBRARIES}")
  include_directories(${SMASH_INCLUDE_DIR})
  # SYSTEM supresses compiling warnings from external libraries
  target_link_libraries(example SYSTEM ${SMASH_LIBRARIES})
  target_link_libraries(example_rate_equations SYSTEM ${SMASH_LIBRARIES})
endif(${SMASH_FOUND})
