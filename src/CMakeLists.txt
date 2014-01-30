include(CheckLibraryExists)

# list the source files
set(smash_src box.cc clebschgordan.cc collisions.cc collider.cc
        cross-sections.cc
        decays.cc distributions.cc experiment.cc fourvector.cc
        input-decaymodes.cc
        input-particles.cc modus.cc param-reader.cc parametrizations.cc
        particles.cc
        outputroutines.cc resonances.cc time.cc)

add_executable(smash smash.cc ${smash_src})
 
# configure a header file to pass some of the CMake settings
# to the source code
# FIXME: hardcoded path
configure_file (
  "${PROJECT_SOURCE_DIR}/src/include/Config.h.in"
  "${PROJECT_BINARY_DIR}/src/include/Config.h"
)
include_directories("${PROJECT_BINARY_DIR}/src")

# directly copy config files to build directory
file(COPY config_general.txt DESTINATION ..)
file(COPY particles.txt DESTINATION ..)
file(COPY decaymodes.txt DESTINATION ..)

# link with gsl and maybe rt
CHECK_LIBRARY_EXISTS(rt clock_gettime "time.h" HAVE_CLOCK_GETTIME)
if(HAVE_CLOCK_GETTIME)
  target_link_libraries(smash gsl gslcblas rt)
else(HAVE_CLOCK_GETTIME)
  target_link_libraries(smash gsl gslcblas)
endif(HAVE_CLOCK_GETTIME)

# tests:
add_subdirectory(tests)
