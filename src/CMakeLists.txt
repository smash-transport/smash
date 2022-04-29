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
# The variable SMASH_VERSION is already set via the project()
# command, since we specify there a VERSION option. Here check
# if .git information is available and increase verbosity
git_describe(SMASH_VERSION_VERBOSE)
if(SMASH_VERSION_VERBOSE)
  set(SMASH_VERSION "${SMASH_VERSION_VERBOSE}")
else()
  set(SMASH_VERSION "${CMAKE_PROJECT_NAME}-${SMASH_VERSION}")
endif()

# Enable warnings for supported compilers only
# --> Here we know the compilers support the flag, do not test!
if(CMAKE_CXX_COMPILER_ID MATCHES "^((Apple)?Clang|GNU)$")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W -Wall -Wextra -Wmissing-declarations -Wpointer-arith -Wshadow -Wuninitialized -Winit-self -Wundef -Wcast-align -Wformat=2 -Wold-style-cast -Werror=switch")
  if (CMAKE_COMPILER_IS_GNUCC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wlogical-op")
    if (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 6.0)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wnull-dereference -Wduplicated-cond")
    endif()
    if (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 7.0)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wrestrict -Wduplicated-branches")
    endif()
  endif()
  if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-missing-braces")
  endif()
  # Enable warnings about deprecated C functions
  option(DEPRECATE_C_FNS "Turn this on to enable warning about deprecated C funtions in SMASH." OFF)
  if (DEPRECATE_C_FNS)
    add_definitions(-include smash/deprecate_c_functions.h -DDONT_ASSERT_TEST)
  endif()
endif()

find_package (Git)
if(GIT_FOUND)
  #message("git found: ${GIT_EXECUTABLE} in version ${GIT_VERSION_STRING}")
  get_filename_component(GIT_PARENT_DIR ${CMAKE_CURRENT_SOURCE_DIR} PATH)
  set(GIT_DIR "${GIT_PARENT_DIR}/.git")
  if(EXISTS "${GIT_DIR}")
    execute_process(
      COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD
      WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
      OUTPUT_VARIABLE GIT_BRANCH)
    string(STRIP ${GIT_BRANCH} GIT_BRANCH)
  else()
    set(GIT_BRANCH "NOTFOUND")
  endif()
  #message("git branch: ${GIT_BRANCH}")
endif (GIT_FOUND)

string(TIMESTAMP BUILD_DATE)

find_package(GSL 2.0 REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(Boost 1.49.0 REQUIRED COMPONENTS filesystem system)

option(TRY_USE_ROOT "Turn this off to disable ROOT output support in SMASH." ON)
if(TRY_USE_ROOT)
  find_package(ROOT 5.34 QUIET)
  if(ROOT_FOUND)
    message(STATUS "Found ROOT ${ROOT_VERSION} (include at ${ROOT_INCLUDE_DIR}).")
    include_directories(SYSTEM "${ROOT_INCLUDE_DIR}")
    set(SMASH_LIBRARIES ${ROOT_LIBRARIES})
    add_definitions(-DSMASH_USE_ROOT)
  else()
    message(STATUS "ROOT not found. Support disabled. See README for more information.")
  endif()
endif()

option(TRY_USE_HEPMC "Turn this off to disable HepMC output support in SMASH." ON)
if(TRY_USE_HEPMC OR TRY_USE_RIVET)
  find_package(HepMC3 QUIET)
  if(HepMC3_FOUND)
    message(STATUS "Found HepMC ${HepMC3_VERSION} (include at ${HEPMC3_INCLUDE_DIR}).")
    include_directories(SYSTEM "${HEPMC3_INCLUDE_DIR}")
    set(SMASH_LIBRARIES
        ${SMASH_LIBRARIES}
        ${HEPMC3_LIBRARIES}
    )
    add_definitions(-DSMASH_USE_HEPMC)
    if((NOT HEPMC3_ROOTIO_LIB) OR (NOT ROOT_FOUND))
        message(STATUS "HepMC3_treeroot output disabled, missing ROOT or/and "
            "HepMC3 ROOT Io support")
    else()
        message(STATUS "HepMC3_treeroot output support enabled")
        add_definitions(-DSMASH_USE_HEPMC_ROOTIO)
    endif()
  else()
    message(STATUS "HepMC not found. Support disabled. See README for more information.")
  endif()
endif()

option(TRY_USE_RIVET "Turn this off to disable Rivet output support in SMASH." ON)
if(TRY_USE_RIVET AND HepMC3_FOUND)
  find_package(Rivet QUIET)
  if(Rivet_FOUND)
    message(STATUS "Found Rivet ${RIVET_VERSION} (include at ${RIVET_INCLUDE_DIR}).")
    include_directories(SYSTEM "${RIVET_INCLUDE_DIR}")
    set(SMASH_LIBRARIES
      ${SMASH_LIBRARIES}
      ${RIVET_LIBRARIES}
      )
    add_definitions(-DSMASH_USE_RIVET)
    set(Rivet_FOUND TRUE)
  else()
    message(STATUS "Rivet not found. Support disabled. See README for more information.")
  endif()
endif()

# find Pythia
find_package(Pythia 8.307 EXACT REQUIRED)
if(Pythia_FOUND)
# the BEFORE at the line below is due to
# https://github.com/smash-transport/smash/issues/38#issuecomment-1033239643
  include_directories(BEFORE SYSTEM
     "${Pythia_INCLUDE_DIRS}"
     )
  set(SMASH_LIBRARIES
     ${SMASH_LIBRARIES}
     ${Pythia_LIBRARIES} ${CMAKE_DL_LIBS}
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
        boxmodus.cc
        binaryoutput.cc
        bremsstrahlungaction.cc
        chemicalpotential.cc
        clebschgordan.cc
        collidermodus.cc
        configuration.cc
        crosssections.cc
        crosssectionsphoton.cc
        customnucleus.cc
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
        fields.cc
        file.cc
        filelock.cc
        fourvector.cc
        fpenvironment.cc
        grandcan_thermalizer.cc
        grid.cc
        hadgas_eos.cc
        hypersurfacecrossingaction.cc
        icoutput.cc
        inputfunctions.cc
        interpolation.cc
        interpolation2D.cc
        isoparticletype.cc
        library.cc
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
        potential_globals.cc
        processbranch.cc
        stringprocess.cc
        propagation.cc
        quantumnumbers.cc
        quantumsampling.cc
        random.cc
        scatteraction.cc
        scatteractionmulti.cc
        scatteractionphoton.cc
        scatteractionsfinder.cc
        setup_particles_decaymodes.cc
        sha256.cc
        spheremodus.cc
        stringfunctions.cc
        tabulation.cc
        thermalizationaction.cc
        thermodynamiclatticeoutput.cc
        thermodynamicoutput.cc
        threevector.cc
        vtkoutput.cc
        wallcrossingaction.cc
        )

if(TRY_USE_ROOT AND ROOT_FOUND)
  set(smash_src ${smash_src} rootoutput.cc)
endif()

if(TRY_USE_HEPMC AND HepMC3_FOUND)
  set(smash_src ${smash_src} hepmcoutput.cc hepmcinterface.cc)
  if(TRY_USE_RIVET AND Rivet_FOUND)
    set(smash_src ${smash_src} rivetoutput.cc)
  endif()
endif()

option(ENABLE_NANOBENCHMARKING "Turn this on to enable code to perform nanobenchmarking in SMASH." OFF)
if(ENABLE_NANOBENCHMARKING)
  if(CMAKE_HOST_SYSTEM_PROCESSOR MATCHES "^arm")
    message(WARNING "Nanobenchmarking is not supported on ARM processors and it was disabled.")
  else()
    set(smash_src ${smash_src} tsc.cc)
    add_definitions(-DUSE_NANOBENCHMARKING_CODE)
  endif()
endif()

# this is the "object library" target: compiles the sources
# only once see https://stackoverflow.com/a/29824424
# NOTE: shared libraries need PIC and this is already set for
#       all targets through CMAKE_POSITION_INDEPENDENT_CODE
add_library(objlib OBJECT ${smash_src})

add_executable(smash smash.cc $<TARGET_OBJECTS:objlib>)

# configure a header file to pass some of the CMake settings
# to the source code
configure_file (
  "${CMAKE_CURRENT_SOURCE_DIR}/include/smash/config.h.in"
  "${CMAKE_CURRENT_BINARY_DIR}/include/smash/config.h"
)

# copy config file to build directory
configure_file(${PROJECT_SOURCE_DIR}/input/config.yaml
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
         -D INPUT_FILE=${PROJECT_SOURCE_DIR}/input/${filename}
         -D OUTPUT_FILE=${CMAKE_CURRENT_BINARY_DIR}/include/${filename}.h
         -D NAME=data
         -P ${PROJECT_SOURCE_DIR}/cmake/copy_file_contents_to_string.cmake
         DEPENDS ${PROJECT_SOURCE_DIR}/input/${filename}
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
   set(CMAKE_REQUIRED_FLAGS ${SANITIZER_FLAG}) # add_compiler_flag only adds the
                                               # flag for compilation - it's
                                               # also required with the link
                                               # command, though
   add_compiler_flag(${SANITIZER_FLAG} CXX_RESULT HAVE_SANITIZER)
   set(CMAKE_REQUIRED_FLAGS "${_tmp}")
   if(NOT HAVE_SANITIZER)
      set(SANITIZER_FLAG -fsanitize=address)
      set(CMAKE_REQUIRED_FLAGS ${SANITIZER_FLAG})
      add_compiler_flag(${SANITIZER_FLAG} CXX_RESULT HAVE_SANITIZER)
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

install(TARGETS smash_shared LIBRARY DESTINATION lib)
install(TARGETS smash RUNTIME DESTINATION bin)

install(DIRECTORY include/smash DESTINATION include)
install(FILES ${generated_headers} DESTINATION include/smash)
install(DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/include/smash" DESTINATION include)
