########################################################
#
#    Copyright (c) 2018-2020,2022
#      SMASH Team
#
#    BSD 3-clause license
#
#########################################################

# cmake-format: off
#=======================================================
# - Try to find SMASH installation
#
# Once done this will define:
#
# SMASH_LIBRARIES      smash libraries
# SMASH_INCLUDE_DIR    directory of smash includes
# SMASH_FOUND          true is smash found
#
# The environment variable SMASH_DIR must be set properly
# to succeed, e.g.:
#   export SMASH_DIR=~/Work/SMASH/smash
#
# Optionally, a SMASH_BUILD_DIR environment variable can
# be defined to specify the path to the folder where SMASH
# has been built. If this is not set, it is assumed that
# SMASH has been built in "${SMASH_DIR}/build".
#=======================================================
# cmake-format: on

# At the moment a macro about the system endianness is needed from within SMASH
include(TestBigEndian)
test_big_endian(IS_BIG_ENDIAN)
if(IS_BIG_ENDIAN)
    message(STATUS "Big endian architecture detected.")
    add_definitions("-DBIG_ENDIAN_ARCHITECTURE")
else()
    message(STATUS "Little endian architecture detected.")
    add_definitions("-DLITTLE_ENDIAN_ARCHITECTURE")
endif()

message(STATUS "Looking for SMASH")

if(NOT DEFINED ENV{SMASH_DIR})
    message(FATAL_ERROR " \n" " Environment variable SMASH_DIR containing the path to the SMASH\n"
                        " codebase found undefined but needed by FindSMASH.cmake module.\n")
endif()

if(NOT DEFINED ENV{SMASH_BUILD_DIR})
    set(SMASH_BUILD_DIR "$ENV{SMASH_DIR}/build")
else()
    set(SMASH_BUILD_DIR "$ENV{SMASH_BUILD_DIR}")
endif()

find_package(GSL 2.0 REQUIRED)
find_package(Eigen3 3.0 REQUIRED)
find_package(Pythia 8.307 EXACT REQUIRED)

set(SMASH_INCLUDE_DIR
    $ENV{SMASH_DIR}/3rdparty/Cuba-4.2.1
    $ENV{SMASH_DIR}/3rdparty/einhard
    $ENV{SMASH_DIR}/3rdparty/yaml-cpp-0.7.0/include
    $ENV{SMASH_DIR}/src/include
    ${SMASH_BUILD_DIR}/src/include # For the decaymodes and particles header files
    ${GSL_INCLUDE_DIR}
    ${EIGEN3_INCLUDE_DIR}
    ${Pythia_INCLUDE_DIRS})
message(VERBOSE "SMASH includes found in ${SMASH_INCLUDE_DIR}")

find_library(SMASH_LIBRARY NAMES smash PATHS ${SMASH_BUILD_DIR}/src)
find_library(EINHARD_LIBRARY NAMES einhard PATHS ${SMASH_BUILD_DIR}/3rdparty/einhard)
find_library(CPPYAML_LIBRARY NAMES yaml-cpp PATHS ${SMASH_BUILD_DIR}/3rdparty/yaml-cpp-0.7.0)
find_library(INTEGRATION_LIBRARY NAMES cuhre PATHS ${SMASH_BUILD_DIR}/3rdparty/Cuba-4.2.1/src/cuhre)
set(SMASH_LIBRARIES
    ${GSL_LIBRARY}
    ${GSL_CBLAS_LIBRARY}
    ${Pythia_LIBRARIES}
    ${CMAKE_DL_LIBS}
    ${EINHARD_LIBRARY}
    ${CPPYAML_LIBRARY}
    ${INTEGRATION_LIBRARY}
    ${SMASH_LIBRARY})
message(VERBOSE "SMASH libraries: ${SMASH_LIBRARIES}")

set(SMASH_FOUND FALSE)
if(SMASH_INCLUDE_DIR
   AND SMASH_LIBRARY
   AND EINHARD_LIBRARY
   AND CPPYAML_LIBRARY
   AND INTEGRATION_LIBRARY)
    set(SMASH_FOUND TRUE)
endif()
message(STATUS "SMASH found: ${SMASH_FOUND}")
