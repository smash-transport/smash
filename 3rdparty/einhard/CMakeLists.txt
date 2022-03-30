cmake_minimum_required(VERSION 3.0)
project(Einhard)

# set the relevant compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fPIC")
if(CXX_MARCH_NATIVE_SUPPORTED)
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
endif()

add_library(einhard STATIC einhard.cpp)

install(FILES einhard.hpp DESTINATION include)
install(FILES stacktrace.h DESTINATION include)
