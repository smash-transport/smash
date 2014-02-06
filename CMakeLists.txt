# The name of our project
project(smash) 

# Minimum cmake version this is tested on
cmake_minimum_required(VERSION 2.8)

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

# have tests run
enable_testing()

# have tests compiled for check target
add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND})

# subdirectories where the code is
add_subdirectory(src) 
