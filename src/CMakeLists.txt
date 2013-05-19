# list the source files
add_executable(mash box.cc collisions.cc fourvector.cc initial-conditions.cc
	param-reader.cc particles.cc outputroutines.cc)

# configure a header file to pass some of the CMake settings
# to the source code
# FIXME: hardcoded path
configure_file (
  "${PROJECT_SOURCE_DIR}/src/include/Config.h.in"
  "${PROJECT_BINARY_DIR}/src/include/Config.h"
)
include_directories("${PROJECT_BINARY_DIR}/src")

# directly copy params file to build directory
file(COPY params.txt DESTINATION ..)

# link with gsl
target_link_libraries(mash gsl gslcblas rt)

# tests:
# verify that the binary has a cli help
add_test(mash_help ${PROJECT_BINARY_DIR}/mash -h)
# verify that the binary runs
add_test(mash_run ${PROJECT_BINARY_DIR}/mash )
# verify that default run has no mem leaks
set(CTEST_MEMORYCHECK_COMMAND valgrind )
set(CTEST_MEMORYCHECK_COMMAND_OPTIONS "-v" )
add_test(memcheck ${CTEST_MEMORYCHECK_COMMAND} ${CTEST_MEMORYCHECK_COMMAND_OPTIONS} ${PROJECT_BINARY_DIR}/mash  "-s 500")
