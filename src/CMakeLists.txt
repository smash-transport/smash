include(CheckLibraryExists)
include(AddCompilerFlag)

find_package(GSL REQUIRED)

# set up include paths
include_directories("${CMAKE_CURRENT_BINARY_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}" "${GSL_INCLUDES}")

# set default libraries for linking
set(SMASH_LIBRARIES ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY})

# if clock_gettime is in rt link against it
CHECK_LIBRARY_EXISTS(rt clock_gettime "time.h" HAVE_CLOCK_GETTIME)
if(HAVE_CLOCK_GETTIME)
   list(APPEND SMASH_LIBRARIES rt)
endif()

# list the source files
set(smash_src boxmodus.cc clebschgordan.cc collisions.cc collidermodus.cc
        cross-sections.cc
        decays.cc distributions.cc experiment.cc fourvector.cc
        input-decaymodes.cc
        input-particles.cc modusdefault.cc
        nucleusmodus.cc
        param-reader.cc parametrizations.cc
        particles.cc
        outputroutines.cc resonances.cc time.cc)

add_executable(smash smash.cc ${smash_src})

# configure a header file to pass some of the CMake settings
# to the source code
configure_file (
  "${CMAKE_CURRENT_SOURCE_DIR}/include/config.h.in"
  "${CMAKE_CURRENT_BINARY_DIR}/include/config.h"
)

# directly copy config files to build directory
macro(copy_at_make)
   foreach(filename ${ARGN})
      add_custom_command(OUTPUT ${PROJECT_BINARY_DIR}/${filename}
         COMMAND ${CMAKE_COMMAND} -E copy
         ${CMAKE_CURRENT_SOURCE_DIR}/${filename}
         ${PROJECT_BINARY_DIR}/${filename}
         DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${filename}
         COMMENT "Copy ${filename} to ${PROJECT_BINARY_DIR}"
         VERBATIM)
      add_custom_target(copy_${filename} ALL DEPENDS ${PROJECT_BINARY_DIR}/${filename})
   endforeach()
endmacro()
copy_at_make(config_general.txt particles.txt decaymodes.txt)

target_link_libraries(smash ${SMASH_LIBRARIES})

# tests:
if(BUILD_TESTING)
   set(ASAN_FLAG -fsanitize=address)
   set(_tmp "${CMAKE_REQUIRED_FLAGS}")
   set(CMAKE_REQUIRED_FLAGS ${ASAN_FLAG}) # AddCompilerFlag only adds the
                                                # flag for compilation - it's
                                                # also required with the link
                                                # command, though
   AddCompilerFlag(${ASAN_FLAG} CXX_RESULT HAVE_ADDRESS_SANITIZER)
   set(CMAKE_REQUIRED_FLAGS "${_tmp}")
   if(HAVE_ADDRESS_SANITIZER)
      add_executable(smash_asan smash.cc ${smash_src})
      set_target_properties(smash_asan PROPERTIES
         COMPILE_FLAGS ${ASAN_FLAG}
         LINK_FLAGS ${ASAN_FLAG}
         )
      target_link_libraries(smash_asan ${SMASH_LIBRARIES})
   endif()

   add_subdirectory(tests)
endif()
