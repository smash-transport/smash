# The name of our project
project(box) 

# Minimum cmake version this is tested on
cmake_minimum_required(VERSION 2.8)

# The current version number
set(VERSION_MAJOR 11042013)

# Set the relevant generic g++ flags (optimisation + warnings)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W -Wall -Wextra -Wmissing-declarations -Wpointer-arith -Wshadow -Wuninitialized  -Winline")

# enable debug preprocessor variable
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG")

# subdirectories where the code is
add_subdirectory(src) 
