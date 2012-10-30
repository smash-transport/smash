# The name of our project
project (NUrQMD) 

# Minimum cmake version this is tested on
cmake_minimum_required(VERSION 2.8)

# The current version number
set (VERSION_MAJOR 0)
set (VERSION_MINOR 0)

# Set all available gcc debug rules
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -W -Wall -Wextra -Wmissing-declarations -Wpointer-arith -Wshadow -Wuninitialized  -Winline")

# subdirectories where the code is
add_subdirectory (src) 
