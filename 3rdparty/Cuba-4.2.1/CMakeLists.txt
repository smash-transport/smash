# This project is added via the add_subdirectory() command
# which creates a new scope and the CMAKE_MINIMUM_REQUIRED_VERSION
# contains the most recent minimum required version in the current
# variable scope or any parent variable scope. Hence, using it
# here should imply that the same SMASH requirement is imposed.
# Omitting the command would be equivalent in this case.
cmake_minimum_required(VERSION ${CMAKE_MINIMUM_REQUIRED_VERSION})
project(Cuba)

include(CheckFunctionExists)
include(CheckIncludeFile)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fomit-frame-pointer -Wall -DNOUNDERSCORE -DREALSIZE=8 -fPIC")
# The '-march=native' flag is not supported e.g. on Apple M1 machines
#  -> check if the C compiler supports it
set(CMAKE_REQUIRED_FLAGS "-march=native")
include(CheckCSourceCompiles)
CHECK_C_SOURCE_COMPILES("int main(void) {return 0;}" C_MARCH_NATIVE_SUPPORTED)
if(C_MARCH_NATIVE_SUPPORTED)
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=native")
else()
   message(STATUS "Your C compiler does not support the '-march=native' flag. This will not be used.")
endif()
# It is crucial to unset the tested flag in order to avoid
# affecting the following calls to 'check_function_exists'
# and 'check_include_file' functions
unset(CMAKE_REQUIRED_FLAGS)

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

install(FILES cuba.h DESTINATION include)
