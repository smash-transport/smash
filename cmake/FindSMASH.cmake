#=======================================================
# - Try to find SMASH instalation
#
# Once done this will define:
#
# SMASH_LIBRARIES      smash libraries
# SMASH_INCLUDE_DIR    directory of smash includes
# SMASH_FOUND          true is smash found
#
# The environment variable SMASH_DIR must be set properly to succeed, e. g.:
# export SMASH_DIR=~/Work/SMASH/smash
#=======================================================

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

message(STATUS "Looking for SMASH ...")

find_package(GSL 2.0 REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(Pythia 8.307 EXACT REQUIRED)

set(SMASH_LIBRARIES ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} ${Pythia_LIBRARIES} -ldl)

set(SMASH_INCLUDE_DIR
    $ENV{SMASH_DIR}/3rdparty/Cuba-4.2.1
    $ENV{SMASH_DIR}/3rdparty/einhard
    $ENV{SMASH_DIR}/3rdparty/yaml-cpp-0.7.0/include
    $ENV{SMASH_DIR}/build/src/include
    $ENV{SMASH_DIR}/src/include
    ${GSL_INCLUDE_DIR}
    ${EIGEN3_INCLUDE_DIR}
    ${Pythia_INCLUDE_DIRS})

message(STATUS "SMASH includes found in ${SMASH_INCLUDE_DIR}")
find_library(SMASH_LIBRARY NAMES smash PATHS $ENV{SMASH_DIR}/build/src)
find_library(EINHARD_LIBRARY NAMES einhard PATHS $ENV{SMASH_DIR}/build/3rdparty/einhard)
find_library(CPPYAML_LIBRARY NAMES yaml-cpp PATHS $ENV{SMASH_DIR}/build/3rdparty/yaml-cpp-0.7.0)
find_library(INTEGRATION_LIBRARY NAMES cuhre
             PATHS $ENV{SMASH_DIR}/build/3rdparty/Cuba-4.2.1/src/cuhre)
set(SMASH_LIBRARIES
    ${SMASH_LIBRARIES}
    ${EINHARD_LIBRARY}
    ${CPPYAML_LIBRARY}
    ${INTEGRATION_LIBRARY}
    ${SMASH_LIBRARY})

# message(STATUS "SMASH libraries: ${SMASH_LIBRARIES}")

set(SMASH_FOUND FALSE)
if(SMASH_INCLUDE_DIR
   AND SMASH_LIBRARY
   AND EINHARD_LIBRARY
   AND CPPYAML_LIBRARY
   AND INTEGRATION_LIBRARY)
    set(SMASH_FOUND TRUE)
endif()

message(STATUS "SMASH found: ${SMASH_FOUND}")
