# CMake minimum version inherited from SMASH context
project(Einhard)

add_library(einhard STATIC einhard.cpp)

if(SMASH_ARCH_FLAG)
    target_add_compiler_flag_if_supported(TARGETS einhard
                                          LANGUAGE CXX
                                          FLAGS "-march=${SMASH_ARCH_FLAG}"
                                          SCOPE PRIVATE)
endif()

set_target_properties(einhard PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED ON)

target_sources(einhard
               PUBLIC
               FILE_SET HEADERS
               BASE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}
               FILES einhard.hpp stacktrace.h)

target_include_directories(einhard SYSTEM PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                                                 $<INSTALL_INTERFACE:include/einhard>)
