cmake_minimum_required(VERSION 2.6.2)
project(Einhard)

# set the relevant compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -march=native -mfpmath=sse -fPIC")

add_library(einhard STATIC einhard.cpp)

install(FILES einhard.hpp TYPE INCLUDE)
install(FILES stacktrace.h TYPE INCLUDE)
