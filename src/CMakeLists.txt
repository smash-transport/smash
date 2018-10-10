########################################################
#
#    Copyright (c) 2014-2015
#      SMASH Team
#
#    BSD 3-clause license
#
#########################################################

include(CheckLibraryExists)
include(AddCompilerFlag)
include(FindPkgConfig)

include(GetGitRevisionDescription)
git_describe(VERSION_MAJOR)

# Enable warnings
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W -Wall -Wextra -Wmissing-declarations -Wpointer-arith -Wshadow -Wuninitialized -Winit-self -Wundef -Wcast-align -Wformat=2 -Wold-style-cast -Werror=switch")
if (CMAKE_COMPILER_IS_GNUCC AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 6.0)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wnull-dereference -Wduplicated-cond")
endif()
if (CMAKE_COMPILER_IS_GNUCC AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 7.0)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wrestrict -Wduplicated-branches")
endif()

find_package (Git)
if(GIT_FOUND)
  #message("git found: ${GIT_EXECUTABLE} in version ${GIT_VERSION_STRING}")
  execute_process(
    COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_BRANCH)
  string(STRIP ${GIT_BRANCH} GIT_BRANCH)
  #message("git branch: ${GIT_BRANCH}")
endif (GIT_FOUND)

string(TIMESTAMP BUILD_DATE)

find_package(GSL 1.15 REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(Boost 1.49.0 REQUIRED COMPONENTS filesystem system)

option(USE_ROOT "Turn this off to disable ROOT output support in SMASH." ON)
if(USE_ROOT)
  find_package(ROOT 5.34)
  if(ROOT_FOUND)
    include_directories(SYSTEM "${ROOT_INCLUDE_DIR}")
    set(SMASH_LIBRARIES ${ROOT_LIBRARIES})
    add_definitions(-DSMASH_USE_ROOT)
  endif()
endif()

# find Pythia
find_package(Pythia REQUIRED)
if(Pythia_FOUND)
  include_directories(SYSTEM
     "${Pythia_INCLUDE_DIRS}"
  )
  set(SMASH_LIBRARIES
     ${SMASH_LIBRARIES}
     ${Pythia_LIBRARIES} -ldl
  )
  add_definitions(
  "-DPYTHIA_XML_DIR=\"${Pythia_xmldoc_PATH}\""
  -DPYTHIA_FOUND
  )
endif()

# set up include paths
include_directories(include)
include_directories(
    "${CMAKE_CURRENT_BINARY_DIR}/include"
    )
include_directories(SYSTEM
    ${GSL_INCLUDE_DIR}
    ${Boost_INCLUDE_DIRS}
    ${EIGEN3_INCLUDE_DIR}
    )


# set default libraries for linking
set(SMASH_LIBRARIES ${SMASH_LIBRARIES}
   ${GSL_LIBRARY}
   ${GSL_CBLAS_LIBRARY}
   ${Boost_LIBRARIES}
   einhard
   yaml-cpp
   cuhre suave divonne vegas  # Cuba multidimensional integration
   )

# list the source files
set(smash_src
        action.cc
        action_globals.cc
        adaptiverejectionsampler.cc
        boxmodus.cc
        binaryoutputcollisions.cc
        binaryoutputparticles.cc
        clebschgordan.cc
        collidermodus.cc
        configuration.cc
        crosssections.cc
        crosssectionsphoton.cc
        decayaction.cc
        decayactionsfinder.cc
        decaymodes.cc
        decaytype.cc
        deformednucleus.cc
        density.cc
        decayactiondilepton.cc
        decayactionsfinderdilepton.cc
        distributions.cc
        energymomentumtensor.cc
        experiment.cc
        file.cc
        filelock.cc
        fourvector.cc
        fpenvironment.cc
        grandcan_thermalizer.cc
        grid.cc
        hadgas_eos.cc
        inputfunctions.cc
        interpolation.cc
        isoparticletype.cc
        listmodus.cc
        logging.cc
        nucleus.cc
        oscaroutput.cc
        pauliblocking.cc
        parametrizations.cc
        particledata.cc
        particles.cc
        particletype.cc
        pdgcode.cc
        potentials.cc
        processbranch.cc
        processstring.cc
        propagation.cc
        quantumnumbers.cc
        random.cc
        scatteraction.cc
        scatteractionphoton.cc
        scatteractionsfinder.cc
        spheremodus.cc
        stringfunctions.cc
        tabulation.cc
        thermalizationaction.cc
        thermodynamicoutput.cc
        threevector.cc
        tsc.cc
        vtkoutput.cc
        wallcrossingaction.cc
        )

if(USE_ROOT AND ROOT_FOUND)
  set(smash_src ${smash_src} rootoutput.cc)
endif()

# this is the "object library" target: compiles the sources only once
# see stackoverflow.com/questions/2152077/is-it-possible-to-get-cmake-to
# -build-both-a-static-and-shared-version-of-the-sam
add_library(objlib OBJECT ${smash_src})
# shared libraries need PIC
set_property(TARGET objlib PROPERTY POSITION_INDEPENDENT_CODE 1)

add_executable(smash smash.cc $<TARGET_OBJECTS:objlib>)

# configure a header file to pass some of the CMake settings
# to the source code
configure_file (
  "${CMAKE_CURRENT_SOURCE_DIR}/include/smash/config.h.in"
  "${CMAKE_CURRENT_BINARY_DIR}/include/smash/config.h"
)

# copy config file to build directory
configure_file(${CMAKE_SOURCE_DIR}/input/config.yaml
               ${PROJECT_BINARY_DIR}/config.yaml
               COPYONLY)
# generate headers from input files
set(generated_headers)
macro(generate_headers)
   foreach(filename ${ARGN})
      list(APPEND generated_headers "${CMAKE_CURRENT_BINARY_DIR}/include/${filename}.h")
      add_custom_command(
         OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/include/${filename}.h
         COMMAND ${CMAKE_COMMAND}
         -D INPUT_FILE=${CMAKE_SOURCE_DIR}/input/${filename}
         -D OUTPUT_FILE=${CMAKE_CURRENT_BINARY_DIR}/include/${filename}.h
         -D NAME=data
         -P ${PROJECT_SOURCE_DIR}/cmake/copy_file_contents_to_string.cmake
         DEPENDS ${CMAKE_SOURCE_DIR}/input/${filename}
         COMMENT "Convert ${CMAKE_SOURCE_DIR}/input/${filename} to ${CMAKE_CURRENT_BINARY_DIR}/include/${filename}.h"
         VERBATIM)
      add_custom_target(generate_${filename}.h ALL DEPENDS
         ${CMAKE_CURRENT_BINARY_DIR}/include/${filename}.h)
   endforeach()
endmacro()
generate_headers(particles.txt decaymodes.txt)

set_source_files_properties(experiment.cc PROPERTIES OBJECT_DEPENDS "${generated_headers}")

target_link_libraries(smash ${SMASH_LIBRARIES})


# Create a shared library out of the whole SMASH
add_library(smash_shared SHARED $<TARGET_OBJECTS:objlib>)
target_link_libraries(smash_shared ${SMASH_LIBRARIES})
set_target_properties(smash_shared PROPERTIES OUTPUT_NAME smash)

# tests:
if(BUILD_TESTING)
   # library for unit tests
   add_library(smash_static STATIC $<TARGET_OBJECTS:objlib>)
   set_target_properties(smash_static PROPERTIES OUTPUT_NAME smash)
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
   set(CMAKE_REQUIRED_FLAGS "${_tmp}")
   if(NOT HAVE_SANITIZER)
      set(SANITIZER_FLAG -fsanitize=address)
      set(CMAKE_REQUIRED_FLAGS ${SANITIZER_FLAG})
      AddCompilerFlag(${SANITIZER_FLAG} CXX_RESULT HAVE_SANITIZER)
      set(CMAKE_REQUIRED_FLAGS "${_tmp}")
   endif()

   if(HAVE_SANITIZER)
      option(USE_SANITIZER "Use available address and undefined behavior sanitizers for test code." ON)
   else()
      set(USE_SANITIZER FALSE)
   endif()

   # define BUILD_TESTS for smash_static in order to enable small code changes for the unit tests.
   # It is an unfortunate side effect that the sanitzer build now also compiles to slightly different code.
   # To fix this we'd have to compile all source files a third time... not sure we prefer that.
   if(USE_SANITIZER)
      set_target_properties(smash_static PROPERTIES
         COMPILE_FLAGS "${SANITIZER_FLAG} -DBUILD_TESTS"
         LINK_FLAGS ${SANITIZER_FLAG}
         )

      add_executable(smash_sanitizer smash.cc)
      set_target_properties(smash_sanitizer PROPERTIES
         COMPILE_FLAGS "${SANITIZER_FLAG} -DBUILD_TESTS"
         LINK_FLAGS ${SANITIZER_FLAG}
         )
      target_link_libraries(smash_sanitizer smash_static ${SMASH_LIBRARIES})
   else()
      set_target_properties(smash_static PROPERTIES COMPILE_FLAGS "-DBUILD_TESTS")
   endif()

   add_subdirectory(tests)
endif()
