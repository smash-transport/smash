include(CheckLibraryExists)
include(AddCompilerFlag)
include(FindPkgConfig)

find_package(GSL REQUIRED)
find_package(Boost 1.49.0 REQUIRED COMPONENTS filesystem system)
pkg_check_modules(YAML REQUIRED yaml-cpp>=0.5.1)

# set up include paths
include_directories(
    "${CMAKE_CURRENT_BINARY_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}"
    )
include_directories(SYSTEM
    "${GSL_INCLUDES}"
    ${Boost_INCLUDE_DIRS}
    ${YAML_INCLUDE_DIRS}
    )

# set default libraries for linking
set(SMASH_LIBRARIES
   ${GSL_LIBRARY}
   ${GSL_CBLAS_LIBRARY}
   ${Boost_LIBRARIES}
   ${YAML_LDFLAGS}
   )

# if clock_gettime is in rt link against it
CHECK_LIBRARY_EXISTS(rt clock_gettime "time.h" HAVE_CLOCK_GETTIME)
if(HAVE_CLOCK_GETTIME)
   list(APPEND SMASH_LIBRARIES rt)
endif()

# list the source files
set(smash_src action.cc boxmodus.cc
        clebschgordan.cc collidermodus.cc configuration.cc cross-sections.cc
        decayaction.cc decayactionsfinder.cc distributions.cc
        experiment.cc fourvector.cc
        modusdefault.cc
        nucleus.cc
        nucleusmodus.cc
        outputroutines.cc
        parametrizations.cc
        particles.cc
        particlesoutput.cc
        random.cc
        resonances.cc
        scatteraction.cc scatteractionsfinder.cc
        time.cc
        vtkoutput.cc
        )

add_executable(smash smash.cc ${smash_src})

# configure a header file to pass some of the CMake settings
# to the source code
configure_file (
  "${CMAKE_CURRENT_SOURCE_DIR}/include/config.h.in"
  "${CMAKE_CURRENT_BINARY_DIR}/include/config.h"
)

# directly copy config files to build directory
configure_file(config.yaml ${PROJECT_BINARY_DIR}/config.yaml COPYONLY)
set(generated_headers)
macro(copy_at_make)
   foreach(filename ${ARGN})
      configure_file(
         ${CMAKE_CURRENT_SOURCE_DIR}/${filename}
         ${PROJECT_BINARY_DIR}/${filename}
         COPYONLY)
      list(APPEND generated_headers "${CMAKE_CURRENT_BINARY_DIR}/${filename}.h")
      add_custom_command(
         OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${filename}.h
         COMMAND ${CMAKE_COMMAND}
         -D INPUT_FILE=${CMAKE_CURRENT_SOURCE_DIR}/${filename}
         -D OUTPUT_FILE=${CMAKE_CURRENT_BINARY_DIR}/${filename}.h
         -D NAME=data
         -P ${PROJECT_SOURCE_DIR}/cmake/copy_file_contents_to_string.cmake
         DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${filename}
         COMMENT "Convert ${CMAKE_CURRENT_SOURCE_DIR}/${filename} to ${CMAKE_CURRENT_BINARY_DIR}/${filename}.h"
         VERBATIM)
      add_custom_target(copy_${filename} ALL DEPENDS
         ${CMAKE_CURRENT_BINARY_DIR}/${filename}.h)
   endforeach()
endmacro()
copy_at_make(particles.txt decaymodes.txt)

set_source_files_properties(experiment.cc PROPERTIES OBJECT_DEPENDS "${generated_headers}")

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
      add_library(smash_lib STATIC ${smash_src})
      set_target_properties(smash_lib PROPERTIES
         COMPILE_FLAGS ${ASAN_FLAG}
         LINK_FLAGS ${ASAN_FLAG}
         )

      add_executable(smash_asan smash.cc)
      set_target_properties(smash_asan PROPERTIES
         COMPILE_FLAGS ${ASAN_FLAG}
         LINK_FLAGS ${ASAN_FLAG}
         )
      target_link_libraries(smash_asan smash_lib ${SMASH_LIBRARIES})

   endif()

   add_subdirectory(tests)
endif()
