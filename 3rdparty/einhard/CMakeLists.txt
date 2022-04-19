# CMake minimum version inherited from SMASH context
project(Einhard)

add_compiler_flags_if_supported(CXX_FLAGS CMAKE_CXX_FLAGS "-std=c++11" "-march=native")

add_library(einhard STATIC einhard.cpp)

install(FILES einhard.hpp DESTINATION include)
install(FILES stacktrace.h DESTINATION include)
