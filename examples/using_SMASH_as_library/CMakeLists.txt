cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(SMASH_as_a_library_example)

# Tell cmake where to find modules
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

add_executable(example example.cc)
include_directories(include)

# Set the relevant generic compiler flags (optimisation + warnings)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W -Wall -Wextra -Wmissing-declarations -Wpointer-arith -Wshadow -Wuninitialized -Winit-self -Wundef -Wcast-align -Wformat=2 -Wold-style-cast -Werror=switch -std=c++11 -mfpmath=sse")
message(STATUS "CXX_FLAGS = " ${CMAKE_CXX_FLAGS})


find_package(SMASH)
if(${SMASH_FOUND})
  include_directories(${SMASH_INCLUDE_DIR})
  target_link_libraries(example ${SMASH_LIBRARIES})
endif(${SMASH_FOUND})
