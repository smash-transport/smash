# list the source files
add_executable(mash box.cc collisions.cc decays.cc distributions.cc
	fourvector.cc initial-conditions.cc input-particles.cc
	param-reader.cc particles.cc propagation.cc outputroutines.cc)

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
file(COPY particles.txt DESTINATION ..)

# link with gsl
target_link_libraries(mash gsl gslcblas rt)

# tests:
add_subdirectory(tests)
