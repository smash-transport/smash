# CMake minimum version inherited from SMASH context

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
   message(FATAL_ERROR "You don't want to configure in the source directory!")
endif()

project(virtest)
install(DIRECTORY vir DESTINATION include/vir)

enable_testing()
add_subdirectory(tests)
