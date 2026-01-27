# CMake minimum version inherited from SMASH context
project(Einhard)

add_library(einhard STATIC einhard.cpp)

target_add_compiler_flag_if_supported(
  TARGETS   einhard
  LANGUAGE  CXX
  FLAGS     "-std=c++11" "-march=native"
  SCOPE     PRIVATE
)

target_sources(einhard PUBLIC
    FILE_SET HEADERS
    BASE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}
    FILES einhard.hpp stacktrace.h
)

target_include_directories(einhard SYSTEM
    PUBLIC
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
      $<INSTALL_INTERFACE:include/einhard>
)
