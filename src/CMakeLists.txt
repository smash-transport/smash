include(CheckLibraryExists)
include(AddCompilerFlag)
include(FindPkgConfig)

include(GetGitRevisionDescription)
git_describe(VERSION_MAJOR)

find_package(GSL REQUIRED)
find_package(Boost 1.49.0 REQUIRED COMPONENTS filesystem system)
pkg_check_modules(YAML REQUIRED yaml-cpp>=0.5.1)

option(USE_ROOT "Turn this on to build ROOT output support into SMASH." OFF)
if(USE_ROOT)
   find_package(ROOT 5.34 REQUIRED)
   include_directories(SYSTEM
      "${ROOT_INCLUDE_DIR}"
      )
   set(SMASH_LIBRARIES
      ${ROOT_LIBRARIES}
      )
   add_definitions(-DSMASH_USE_ROOT)
endif()

# set up include paths
include_directories(
    "${CMAKE_CURRENT_BINARY_DIR}"
    )
include_directories(SYSTEM
    "${GSL_INCLUDES}"
    ${Boost_INCLUDE_DIRS}
    ${YAML_INCLUDE_DIRS}
    )

# set default libraries for linking
set(SMASH_LIBRARIES ${SMASH_LIBRARIES}
   ${GSL_LIBRARY}
   ${GSL_CBLAS_LIBRARY}
   ${Boost_LIBRARIES}
   ${YAML_LDFLAGS}
   einhard
   )

# list the source files
set(smash_src
        action.cc
        boxmodus.cc
        binaryoutputcollisions.cc
        binaryoutputparticles.cc
        collidermodus.cc
        configuration.cc
        decayaction.cc
        decayactionsfinder.cc
        decaymodes.cc
        deformednucleus.cc
        distributions.cc
        experiment.cc
        fourvector.cc
        inputfunctions.cc
        logging.cc
        modusdefault.cc
        nucleus.cc
        nucleusmodus.cc
        oscaroutput.cc
        parametrizations.cc
        particledata.cc
        particles.cc
        particletype.cc
        pdgcode.cc
        processbranch.cc
        random.cc
        resonances.cc
        scatteraction.cc
        scatteractionsfinder.cc
        spheremodus.cc
        stringfunctions.cc
        threevector.cc
        vtkoutput.cc
        width.cc
        )

if(USE_ROOT)
   set(smash_src ${smash_src}
      rootoutput.cc
      )
endif()

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
   # library for unit tests
   add_library(smash_lib STATIC ${smash_src})

   # Look for sanitizers supported by the compiler. Sanitizers add
   # instrumentation code to the executable that test correctness of the code:
   # * address: looks for use-after-free and out-of-bounds like errors
   # * undefined: looks for undefined behavior in the code like misaligned
   #              accesses, creation of infinity/NaN values, or conversions to
   #              unrepresentable values
   # If the undefined sanitizer is not available with the compiler it will fall
   # back to address sanitizer only. If this ones also not available the
   # sanitizer build will be skipped altogether.
   set(SANITIZER_FLAG -fsanitize=address,undefined)
   set(_tmp "${CMAKE_REQUIRED_FLAGS}")
   set(CMAKE_REQUIRED_FLAGS ${SANITIZER_FLAG}) # AddCompilerFlag only adds the
                                               # flag for compilation - it's
                                               # also required with the link
                                               # command, though
   AddCompilerFlag(${SANITIZER_FLAG} CXX_RESULT HAVE_SANITIZER)
   if(NOT HAVE_SANITIZER)
      set(SANITIZER_FLAG -fsanitize=address)
      AddCompilerFlag(${SANITIZER_FLAG} CXX_RESULT HAVE_SANITIZER)
   endif()
   set(CMAKE_REQUIRED_FLAGS "${_tmp}")
   if(HAVE_SANITIZER)
      set_target_properties(smash_lib PROPERTIES
         COMPILE_FLAGS ${SANITIZER_FLAG}
         LINK_FLAGS ${SANITIZER_FLAG}
         )

      add_executable(smash_sanitizer smash.cc)
      set_target_properties(smash_sanitizer PROPERTIES
         COMPILE_FLAGS ${SANITIZER_FLAG}
         LINK_FLAGS ${SANITIZER_FLAG}
         )
      target_link_libraries(smash_sanitizer smash_lib ${SMASH_LIBRARIES})
   endif()

   add_subdirectory(tests)
endif()
