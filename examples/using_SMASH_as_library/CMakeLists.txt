cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
project(SMASH_as_a_library_example)

# Tell cmake where to find modules (also append ../../cmake for testing)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake"
     "${CMAKE_CURRENT_LIST_DIR}/../../cmake")

add_executable(example example.cc)
add_executable(example_rate_equations example_rate_equations.cc)
add_executable(example_smash_wrapper example_smash_wrapper.cc)

# Set the relevant generic compiler flags, trying to be close to what SMASH uses
set(CMAKE_CXX_FLAGS
    "${CMAKE_CXX_FLAGS} -W -Wall -Wextra -Wmissing-declarations -Wpointer-arith -Wshadow -Wuninitialized -Winit-self -Wundef -Wcast-align -Wformat=2 -Wold-style-cast -Werror=switch"
)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-math-errno -std=c++17")
message(STATUS "CXX_FLAGS = " ${CMAKE_CXX_FLAGS})

find_package(SMASH)
if(SMASH_FOUND)
    message(STATUS "SMASH libraries = ${SMASH_LIBRARIES}")
    include_directories(SYSTEM ${SMASH_INCLUDE_DIR})
    target_link_libraries(example ${SMASH_LIBRARIES})
    target_link_libraries(example_rate_equations ${SMASH_LIBRARIES})
    target_link_libraries(example_smash_wrapper ${SMASH_LIBRARIES})
    add_definitions("-DSMASH_INPUT_DIR=\"${SMASH_INPUT_FILES_DIR}\"")
else()
    message(FATAL_ERROR "SMASH libraries not found!")
endif()
