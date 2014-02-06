include(CheckLibraryExists)

# find GSL
find_path(GSL_INCLUDES gsl/gsl_sf_bessel.h)
find_library(GSL_LIBRARY gsl)
find_library(GSL_CBLAS_LIBRARY gslcblas)

# set up include paths
include_directories("${CMAKE_CURRENT_BINARY_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}" "${GSL_INCLUDES}")

# list the source files
set(smash_src boxmodus.cc clebschgordan.cc collisions.cc collidermodus.cc
        cross-sections.cc
        decays.cc distributions.cc experiment.cc fourvector.cc
        input-decaymodes.cc
        input-particles.cc modusdefault.cc param-reader.cc parametrizations.cc
        particles.cc
        outputroutines.cc resonances.cc time.cc)

add_executable(smash smash.cc ${smash_src})
 
# configure a header file to pass some of the CMake settings
# to the source code
configure_file (
  "${CMAKE_CURRENT_SOURCE_DIR}/include/Config.h.in"
  "${CMAKE_CURRENT_BINARY_DIR}/include/Config.h"
)

# directly copy config files to build directory
file(COPY config_general.txt DESTINATION ..)
file(COPY particles.txt DESTINATION ..)
file(COPY decaymodes.txt DESTINATION ..)

# link with gsl and maybe rt
CHECK_LIBRARY_EXISTS(rt clock_gettime "time.h" HAVE_CLOCK_GETTIME)
if(HAVE_CLOCK_GETTIME)
  target_link_libraries(smash ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} rt)
else(HAVE_CLOCK_GETTIME)
  target_link_libraries(smash ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY})
endif(HAVE_CLOCK_GETTIME)

# tests:
include (CTest)
if(BUILD_TESTING)
   add_subdirectory(tests)
endif()
