cmake_minimum_required(VERSION 2.6.2)
project(Cuba)

include(CheckFunctionExists)
include(CheckIncludeFile)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fomit-frame-pointer -Wall -DNOUNDERSCORE -DREALSIZE=8 -march=native -mfpmath=sse -fPIC")

set(CMAKE_REQUIRED_INCLUDES "unistd.h" "sys/ipc.h" "sys/shm.h" "math.h" "stdlib.h")

set(COMMON_DEPS
    "${CMAKE_CURRENT_SOURCE_DIR}/src/common/Data.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/common/Global.c"
)

check_function_exists(fork HAVE_FORK)
check_function_exists(shmget HAVE_SHMGET)
check_function_exists(erf HAVE_ERF)
check_function_exists(getloadavg HAVE_GETLOADAVG)
check_include_file(alloca.h HAVE_ALLOCA_H)

# Enabling support for parallelization causes compilation to fail.
# Disabled for now.
#if(${HAVE_FORK})
#    add_definitions(-DHAVE_FORK)
#endif()
if(${HAVE_SHMGET})
    add_definitions(-DHAVE_SHMGET)
endif()
if(${HAVE_ERF})
    add_definitions(-DHAVE_ERF)
endif()
if(${HAVE_GETLOADAVG})
    add_definitions(-DHAVE_GETLOADAVG)
endif()
if(${HAVE_ALLOCA_H})
    add_definitions(-DHAVE_ALLOCA_H)
endif()

add_subdirectory(src/cuhre)
add_subdirectory(src/suave)
add_subdirectory(src/vegas)
add_subdirectory(src/divonne)
