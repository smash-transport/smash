########################################################
#
#    Copyright (c) 2017-2018,2021-2023,2026
#      SMASH Team
#
#    BSD 3-clause license
#
#########################################################

# CMake minimum version inherited from SMASH context
project(Cuba LANGUAGES C)

include(CheckFunctionExists)
include(CheckIncludeFile)

set(COMMON_DEPS "${CMAKE_CURRENT_SOURCE_DIR}/src/common/Data.c"
                "${CMAKE_CURRENT_SOURCE_DIR}/src/common/Global.c")

set(CMAKE_REQUIRED_INCLUDES
    "unistd.h"
    "sys/ipc.h"
    "sys/shm.h"
    "math.h"
    "stdlib.h")

check_function_exists(shmget HAVE_SHMGET)
check_function_exists(erf HAVE_ERF)
check_function_exists(getloadavg HAVE_GETLOADAVG)
check_include_file(alloca.h HAVE_ALLOCA_H)

set(CUBA_COMPILE_DEFINITIONS NOUNDERSCORE REALSIZE=8)

# Enabling support for parallelization causes compilation to fail. Disabled for now.
# ~~~
#  check_function_exists(fork HAVE_FORK)
#  if(${HAVE_FORK})
#      list(APPEND CUBA_COMPILE_DEFINITIONS HAVE_FORK)
#  endif()
# ~~~
foreach(_definition "HAVE_SHMGET" "HAVE_ERF" "HAVE_GETLOADAVG" "HAVE_ALLOCA_H")
    if(${_definition})
        list(APPEND CUBA_COMPILE_DEFINITIONS ${_definition})
    endif()
endforeach()

add_subdirectory(src/cuhre)
add_subdirectory(src/suave)
add_subdirectory(src/vegas)
add_subdirectory(src/divonne)

set(CUBA_TARGETS cuhre suave vegas divonne)

target_add_compiler_flag_if_supported(TARGETS ${CUBA_TARGETS}
                                      LANGUAGE C
                                      FLAGS "-Wall" "-march=native"
                                      SCOPE PRIVATE)

foreach(_target IN LISTS CUBA_TARGETS)

    target_include_directories(${_target} SYSTEM
                               PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/${_target}>
                                       $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/common>
                               PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                                      $<INSTALL_INTERFACE:include/cuba>)

    target_sources(${_target}
                   PUBLIC
                   FILE_SET HEADERS
                   BASE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}
                   FILES cuba.h)

    target_compile_definitions(${_target} PRIVATE ${CUBA_COMPILE_DEFINITIONS})

    # Flags -pg and -fomit-frame-pointer are incompatible, hence omit the latter in profiling mode,
    # where performance is irrelevant
    if(NOT CMAKE_BUILD_TYPE STREQUAL "Profiling")
        target_add_compiler_flag_if_supported(TARGETS ${_target}
                                              LANGUAGE C
                                              FLAGS -fomit-frame-pointer
                                              SCOPE PRIVATE)
    endif()

endforeach()
