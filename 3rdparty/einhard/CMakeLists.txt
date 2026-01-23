# CMake minimum version inherited from SMASH context
project(Einhard)

add_library(einhard STATIC einhard.cpp)

target_add_compiler_flag_if_supported(
  TARGETS   einhard
  LANGUAGE  CXX
  FLAGS     "-std=c++11" "-march=native"
  SCOPE     PRIVATE
)
