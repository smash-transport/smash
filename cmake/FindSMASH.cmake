# Try to find SMASH instalation
#
# Once done this will define:
#
# SMASH_LIBRARIES      smash libraries
# SMASH_INCLUDE_DIR    directory of smash includes
# SMASH_FOUND          true is smash found
#
# The environment variable SMASH_DIR must be set properly to succeed, e. g.:
# export SMASH_DIR=~/Work/SMASH/smash

message(STATUS "Looking for SMASH ...")

find_package(GSL 2.0 REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(Boost 1.49.0 REQUIRED COMPONENTS filesystem system)
find_package(Pythia 8.235 EXACT REQUIRED)

set(SMASH_LIBRARIES ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} ${Boost_LIBRARIES} ${Pythia_LIBRARIES} -ldl)

set(SMASH_INCLUDE_DIR
   $ENV{SMASH_DIR}/3rdparty/Cuba-4.2
   $ENV{SMASH_DIR}/3rdparty/einhard
   $ENV{SMASH_DIR}/3rdparty/yaml-cpp-0.6.2/include
   $ENV{SMASH_DIR}/build/src/include
   $ENV{SMASH_DIR}/src/include
   ${GSL_INCLUDE_DIR}
   ${Boost_INCLUDE_DIRS}
   ${EIGEN3_INCLUDE_DIR}
   ${Pythia_INCLUDE_DIRS}
)

message(STATUS "SMASH includes found in ${SMASH_INCLUDE_DIR}")
find_library(SMASH_LIBRARY        NAMES smash     PATHS $ENV{SMASH_DIR}/build/src)
find_library(EINHARD_LIBRARY      NAMES einhard   PATHS $ENV{SMASH_DIR}/build/3rdparty/einhard)
find_library(CPPYAML_LIBRARY      NAMES yaml-cpp  PATHS $ENV{SMASH_DIR}/build/3rdparty/yaml-cpp-0.6.2)
find_library(INTEGRATION_LIBRARY  NAMES cuhre     PATHS $ENV{SMASH_DIR}/build/3rdparty/Cuba-4.2/src/cuhre)
set(SMASH_LIBRARIES ${SMASH_LIBRARIES}
   ${EINHARD_LIBRARY}
   ${CPPYAML_LIBRARY}
   ${INTEGRATION_LIBRARY}
   ${SMASH_LIBRARY}
)

#message(STATUS "SMASH libraries: ${SMASH_LIBRARIES}")

set(SMASH_FOUND FALSE)
if (SMASH_INCLUDE_DIR AND
    SMASH_LIBRARY     AND
    EINHARD_LIBRARY   AND
    CPPYAML_LIBRARY   AND
    INTEGRATION_LIBRARY)
   set(SMASH_FOUND TRUE)
endif ()

message(STATUS "SMASH found: ${SMASH_FOUND}")
