########################################################
#
#    Copyright (c) 2012-2025
#      SMASH Team
#
#    BSD 3-clause license
#
########################################################

include(CheckLibraryExists)
include(AddCompilerFlag)
include(FindPkgConfig)

# According to modern CMake mental model: Declare the targets early, and then feed it sources, link
# dependencies, compile definitions, etc. incrementally. The "object library" target is used to
# compile the sources only once see https://stackoverflow.com/a/29824424 and all its features will
# be PUBLIC since this is a library to serve internally both the executable and the shared library
# and all its properties must propagate. Here we use the smash_api interface library to collect
# everything that should be exposed to users of the shared library. NOTE: shared libraries need PIC
# and this is already set for all targets through CMAKE_POSITION_INDEPENDENT_CODE
add_library(smash_objlib OBJECT)
add_executable(smash smash.cc $<TARGET_OBJECTS:smash_objlib>)
add_library(smash_shared SHARED $<TARGET_OBJECTS:smash_objlib>)
add_library(smash_api INTERFACE)

# Both smash executable and shared library are linked privately to the object library because the
# object is only an implementation detail and we do not want to expose it to clients, but also to
# avoid that the dependency on the external libraries is propagated to users of the shared library,
# which should not be aware of them (dependency on external library is introduced in the Config
# file). What has to be propagated to clients belongs to smash_api interface target.
target_link_libraries(smash_objlib PRIVATE smash_api)
target_link_libraries(smash PRIVATE smash_objlib smash_api)
target_link_libraries(smash_shared PRIVATE smash_objlib PUBLIC smash_api)

# Call shared library smash (e.g. libsmash.so). Note that, since "the build-tree runtime path will
# not be used for binaries in the install tree" (CMake documentation), we request here to set the
# INSTALL_RPATH in the library to ensure that projects linking to SMASH finds at run time libraries
# like Pythia, which are likely to be in not standard places.
set_target_properties(smash_shared PROPERTIES OUTPUT_NAME smash INSTALL_RPATH "${ORIGIN}"
                                              INSTALL_RPATH_USE_LINK_PATH TRUE)

# Prepare list variables to collect libraries, compiler definitions for internal usage in SMASH,
# compiler definitions for users of SMASH, internal compile options and so on
target_link_libraries(smash_objlib
                      PUBLIC einhard
                             yaml-cpp
                             cuhre
                             suave
                             divonne
                             vegas)

# Enable warnings for supported compilers only --> Here we know the compilers support the flag, do
# not test!
if(CMAKE_CXX_COMPILER_ID MATCHES "^((Apple)?Clang|GNU)$")
    set(CMAKE_CXX_FLAGS
        "${CMAKE_CXX_FLAGS} -W -Wall -Wextra -Wmissing-declarations -Wpointer-arith -Wshadow -Wuninitialized -Wfloat-conversion -Winit-self -Wundef -Wcast-align -Wformat=2 -Wold-style-cast -Werror=switch"
    )
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        set(CMAKE_CXX_FLAGS
            "${CMAKE_CXX_FLAGS} -Wlogical-op -Wnull-dereference -Wduplicated-cond -Wrestrict -Wduplicated-branches -Wsuggest-override"
        )
        # Needed for workaround about C++17 compatibility (bug in some GNU versions)
        add_compile_definitions(GCC_COMPILER)
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        set(CMAKE_CXX_FLAGS
            "${CMAKE_CXX_FLAGS} -Wno-missing-braces -Winconsistent-missing-override")
        if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "14")
            # Use of contracted mathematical expression like FMA (if the hardware supports them) are
            # controlled by a flag whose default unfortunately varies among compilers and versions
            # thereof. Since clang-14 changed the default from "off" to "on", we set here the value
            # to "on" also for older versions. GNU default is "fast", but this seems to be similar
            # to clang's "on" (which fuses in the same statement) than to clang's "fast" (which
            # fuses across statements). Here we simply try to have the same behaviour across the
            # supported compilers.
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ffp-contract=on")
        endif()
    endif()
    # Enable warnings about deprecated C functions
    option(DEPRECATE_C_FNS "Turn this on to enable warning about deprecated C functions in SMASH."
           OFF)
    if(DEPRECATE_C_FNS)
        target_compile_definitions(smash_objlib PUBLIC "DONT_ASSERT_TEST")
        target_compile_options(smash_objlib PUBLIC "-include" "smash/deprecate_c_functions.h")
    endif()
endif()

find_package(Git)
if(GIT_FOUND)
    get_filename_component(GIT_PARENT_DIR ${CMAKE_CURRENT_SOURCE_DIR} PATH)
    set(GIT_DIR "${GIT_PARENT_DIR}/.git")
    if(EXISTS "${GIT_DIR}")
        execute_process(COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD
                        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} OUTPUT_VARIABLE GIT_BRANCH)
        string(STRIP ${GIT_BRANCH} GIT_BRANCH)
    else()
        set(GIT_BRANCH "NOTFOUND")
    endif()
endif(GIT_FOUND)

string(TIMESTAMP BUILD_DATE)

find_package(GSL 2.0 REQUIRED)
if(GSL_FOUND)
    if(TARGET GSL::gsl AND TARGET GSL::gslcblas)
        create_system_wrapper_target(GSL_bundle GSL::gsl GSL::gslcblas)
        target_link_libraries(smash_objlib PUBLIC GSL_bundle)
    else()
        message(FATAL_ERROR " \n" " GSL library found, but its targets are somehow missing.\n"
                            " Ensure to be using the official CMake module to find it.\n")
    endif()
endif()

find_package(Eigen3 3.0 REQUIRED)
if(EIGEN3_FOUND)
    if(TARGET Eigen3::Eigen)
        target_link_libraries(smash_objlib PUBLIC Eigen3::Eigen)
    else()
        message(FATAL_ERROR " \n" " Eigen3 library found, but its target is somehow missing.\n"
                            " Ensure to be using the provided CMake module to find it.\n")
    endif()
endif()

option(TRY_USE_ROOT "Turn this off to disable ROOT output support in SMASH." ON)
if(TRY_USE_ROOT)
    find_package(ROOT CONFIG QUIET COMPONENTS Core RIO Tree)
    set(ROOT_MINIMUM_VERSION "6.20")
    if(ROOT_FOUND)
        message(STATUS "Candidate ROOT version ${ROOT_VERSION} found, "
                       "requested ${ROOT_MINIMUM_VERSION}")
        if(ROOT_VERSION VERSION_LESS ROOT_MINIMUM_VERSION)
            set(ROOT_FOUND FALSE)
        endif()
    endif()
    if(ROOT_FOUND)
        message(STATUS "Found ROOT: ${ROOT_INCLUDE_DIRS} (version ${ROOT_VERSION})")
        create_system_wrapper_target(ROOT_bundle ROOT::Core ROOT::RIO ROOT::Tree)
        target_link_libraries(smash_objlib PUBLIC ROOT_bundle)
        target_compile_definitions(smash_api INTERFACE "SMASH_USE_ROOT")
    else()
        message(STATUS "Suitable ROOT not found. Support disabled. "
                       "See README for more information.")
    endif()
endif()

option(TRY_USE_HEPMC "Turn this off to disable HepMC3 output support in SMASH." ON)
if(TRY_USE_HEPMC OR TRY_USE_RIVET)
    find_package(HepMC3 QUIET CONFIG) # Min version not handled here but after to better inform user
    set(HepMC3_MINIMUM_VERSION "3.2.6")
    if(HepMC3_FOUND)
        message(STATUS "Candidate HepMC3 version ${HepMC3_VERSION} found, "
                       "requested ${HepMC3_MINIMUM_VERSION}")
        if(HepMC3_VERSION VERSION_LESS HepMC3_MINIMUM_VERSION)
            set(HepMC3_FOUND FALSE)
        endif()
    endif()
    if(HepMC3_FOUND)
        message(STATUS "Found HepMC3: ${HEPMC3_INCLUDE_DIR} (version ${HepMC3_VERSION})")
        if(NOT TARGET HepMC3::HepMC3)
            message(FATAL_ERROR "HepMC3 found, but its main target is missing.")
        else()
            set(HEPMC3_TARGETS HepMC3::HepMC3)
        endif()
        if(NOT ROOT_FOUND)
            message(STATUS "HepMC3_treeroot output disabled, missing ROOT support")
        elseif(NOT TARGET HepMC3::rootIO)
            message(STATUS "HepMC3_treeroot output disabled, missing HepMC3 ROOT I/O support")
        else()
            message(STATUS "HepMC3_treeroot output support enabled")
            list(APPEND HEPMC3_TARGETS HepMC3::rootIO)
            target_compile_definitions(smash_api INTERFACE "SMASH_USE_HEPMC_ROOTIO")
        endif()
        create_system_wrapper_target(HepMC3_bundle ${HEPMC3_TARGETS})
        target_link_libraries(smash_objlib PUBLIC HepMC3_bundle)
        target_compile_definitions(smash_api INTERFACE "SMASH_USE_HEPMC")
    else()
        message(STATUS "Suitable HepMC3 not found. Support disabled. "
                       "See README for more information.")
    endif()
endif()

option(TRY_USE_RIVET "Turn this off to disable Rivet output support in SMASH." ON)
if(TRY_USE_RIVET AND HepMC3_FOUND)
    find_package(Rivet QUIET) # Min version not handled here but after to better inform user
    set(Rivet_MINIMUM_VERSION "3.1")
    if(Rivet_FOUND)
        message(STATUS "Candidate Rivet version ${RIVET_VERSION} found, "
                       "requested ${Rivet_MINIMUM_VERSION}")
        if(RIVET_VERSION VERSION_LESS Rivet_MINIMUM_VERSION)
            set(Rivet_FOUND FALSE)
        elseif(RIVET_VERSION VERSION_GREATER_EQUAL "3.2")
            message(ATTENTION "Rivet version larger than 3.1.x found. "
                              "The SMASH interface might miss some functionality.")
        endif()
    endif()
    if(Rivet_FOUND)
        message(STATUS "Found Rivet ${RIVET_INCLUDE_DIR} (version ${RIVET_VERSION})")
        if(NOT TARGET Rivet::Rivet)
            message(FATAL_ERROR " \n" " Rivet library found, but its target is somehow missing.\n"
                                " Ensure to be using the provided CMake module to find it.\n")
        endif()
        target_link_libraries(smash_objlib PUBLIC Rivet::Rivet)
        target_compile_definitions(smash_api INTERFACE "SMASH_USE_RIVET")
    else()
        message(STATUS "Suitable Rivet not found. Support disabled. "
                       "See README for more information.")
    endif()
endif()

find_package(Pythia 8.316 EXACT REQUIRED)
if(Pythia_FOUND)
    if(NOT TARGET Pythia::Pythia)
        message(FATAL_ERROR "Pythia library found, but its target is somehow missing.\n"
                            "Ensure to be using the provided CMake module to find it.\n")
    endif()
    # Here we PREPEND the Pythia target to the list of libraries to ensure that it is linked before
    # all the other and it is meant to avoid possible issues with e.g. ROOT whose include
    # directories might contain a Pythia.h header that might not be the one of the Pythia linked to
    # SMASH. See https://github.com/smash-transport/smash/issues/38#issuecomment-1033239643
    target_link_libraries(smash_objlib PUBLIC Pythia::Pythia)
    target_compile_definitions(smash_objlib PUBLIC "PYTHIA_XML_DIR=\"${Pythia_xmldoc_PATH}\""
                                                   "PYTHIA_FOUND")
endif()

# set up include paths
include_directories(include)
include_directories("${CMAKE_CURRENT_BINARY_DIR}/include")

# ~~~
# For GNU and Clang compilers, link to C++ filesystem library if compiler version requires it
#  -> See "Notes" at https://en.cppreference.com/w/cpp/filesystem
#
# NOTE: As of Xcode 11 and MacOS 10.15, <filesystem> is shipped in the main
#       libc++.dylib and there's no need to link any additional library.
# ~~~
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9.1)
        target_link_libraries(smash_objlib PUBLIC stdc++fs)
    endif()
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9.0)
        target_link_libraries(smash_objlib PUBLIC c++fs)
    endif()
endif()

# list the source files
target_sources(smash_objlib
               PRIVATE action.cc
                       boxmodus.cc
                       binaryoutput.cc
                       bremsstrahlungaction.cc
                       chemicalpotential.cc
                       clebschgordan.cc
                       clebschgordan_lookup.cc
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
                       alphaclusterednucleus.cc
                       density.cc
                       decayactiondilepton.cc
                       decayactionsfinderdilepton.cc
                       dynamicfluidfinder.cc
                       distributions.cc
                       energymomentumtensor.cc
                       experiment.cc
                       fields.cc
                       file.cc
                       filelock.cc
                       fluidizationaction.cc
                       fourvector.cc
                       fpenvironment.cc
                       grandcan_thermalizer.cc
                       grid.cc
                       hadgas_eos.cc
                       hypersurfacecrossingfinder.cc
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
                       stringify.cc
                       tabulation.cc
                       thermalizationaction.cc
                       thermodynamiclatticeoutput.cc
                       thermodynamicoutput.cc
                       threevector.cc
                       vtkoutput.cc
                       wallcrossingaction.cc)

if(TRY_USE_ROOT AND ROOT_FOUND)
    target_sources(smash_objlib PRIVATE rootoutput.cc)
endif()

if(TRY_USE_HEPMC AND HepMC3_FOUND)
    target_sources(smash_objlib PRIVATE hepmcoutput.cc hepmcinterface.cc)
    if(TRY_USE_RIVET AND Rivet_FOUND)
        target_sources(smash_objlib PRIVATE rivetoutput.cc)
    endif()
endif()

option(ENABLE_NANOBENCHMARKING "Turn this on to enable code to perform nanobenchmarking in SMASH."
       OFF)
if(ENABLE_NANOBENCHMARKING)
    if(CMAKE_HOST_SYSTEM_PROCESSOR MATCHES "^arm")
        message(WARNING "Nanobenchmarking is not supported on ARM processors and it was disabled.")
    else()
        target_sources(smash_objlib PRIVATE tsc.cc)
        target_compile_definitions(smash_objlib PUBLIC "USE_NANOBENCHMARKING_CODE")
    endif()
endif()

# configure a header file to pass some of the CMake settings to the source code
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/include/smash/config.h.in"
               "${CMAKE_CURRENT_BINARY_DIR}/include/smash/config.h")

# copy config file to build directory
configure_file(${PROJECT_SOURCE_DIR}/input/config.yaml ${PROJECT_BINARY_DIR}/config.yaml COPYONLY)
# generate headers from input files
set(generated_headers)
macro(generate_headers)
    foreach(filename ${ARGN})
        list(APPEND generated_headers "${CMAKE_CURRENT_BINARY_DIR}/include/${filename}.h")
        add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/include/${filename}.h
                           COMMAND ${CMAKE_COMMAND} -D
                                   INPUT_FILE=${PROJECT_SOURCE_DIR}/input/${filename} -D
                                   OUTPUT_FILE=${CMAKE_CURRENT_BINARY_DIR}/include/${filename}.h -D
                                   NAME=data -P
                                   ${PROJECT_SOURCE_DIR}/cmake/copy_file_contents_to_string.cmake
                           DEPENDS ${PROJECT_SOURCE_DIR}/input/${filename}
                           COMMENT "Convert ${CMAKE_SOURCE_DIR}/input/${filename} to ${CMAKE_CURRENT_BINARY_DIR}/include/${filename}.h"
                           VERBATIM)
        add_custom_target(generate_${filename}.h ALL
                          DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/include/${filename}.h)
    endforeach()
endmacro()
generate_headers(particles.txt decaymodes.txt)
set_source_files_properties(experiment.cc PROPERTIES OBJECT_DEPENDS "${generated_headers}")

# tests:
if(BUILD_TESTING)
    # library for unit tests
    add_library(smash_testlib STATIC $<TARGET_OBJECTS:smash_objlib>)
    target_link_libraries(smash_testlib PUBLIC smash_objlib)
    set_target_properties(smash_testlib PROPERTIES OUTPUT_NAME smash)
    # Look for sanitizers supported by the compiler. Sanitizers add instrumentation code to the
    # executable that test correctness of the code:
    #
    # * address: looks for use-after-free and out-of-bounds like errors
    # * undefined: looks for undefined behavior in the code like misaligned accesses, creation of
    #   infinity/NaN values, or conversions to unrepresentable values
    #
    # If the undefined sanitizer is not available with the compiler it will fall back to address
    # sanitizer only. If this ones also not available the sanitizer build will be skipped
    # altogether.
    set(SANITIZER_FLAG -fsanitize=address,undefined)
    set(_tmp "${CMAKE_REQUIRED_FLAGS}")
    set(CMAKE_REQUIRED_FLAGS ${SANITIZER_FLAG}) # add_compiler_flag only adds the flag for
                                                # compilation - it's also required with the link
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
        option(USE_SANITIZER
               "Use available address and undefined behavior sanitizers for test code." ON)
    else()
        set(USE_SANITIZER FALSE)
    endif()

    # define BUILD_TESTS for smash_testlib in order to enable small code changes for the unit tests.
    # It is an unfortunate side effect that the sanitzer build now also compiles to slightly
    # different code. To fix this we'd have to compile all source files a third time... not sure we
    # prefer that.
    if(USE_SANITIZER)
        set_target_properties(smash_testlib
                              PROPERTIES COMPILE_FLAGS "${SANITIZER_FLAG} -DBUILD_TESTS"
                                         LINK_FLAGS ${SANITIZER_FLAG})

        add_executable(smash_sanitizer smash.cc)
        set_target_properties(smash_sanitizer
                              PROPERTIES COMPILE_FLAGS "${SANITIZER_FLAG} -DBUILD_TESTS"
                                         LINK_FLAGS ${SANITIZER_FLAG})
        target_link_libraries(smash_sanitizer smash_testlib)
    else()
        set_target_properties(smash_testlib PROPERTIES COMPILE_FLAGS "-DBUILD_TESTS")
    endif()

    add_subdirectory(tests)
endif()

# Since at the top level the dependency of the install target on all target was disabled, it is
# needed to now build the smash library and exectuable target here before installing them
install(CODE "execute_process(COMMAND ${CMAKE_COMMAND}
                                      --build ${PROJECT_BINARY_DIR}
                                      --target smash_shared smash)")
install(TARGETS smash_shared LIBRARY DESTINATION "lib")
install(FILES ${generated_headers} DESTINATION "include/smash")
# Note that the absent trailing slash in the DIRECTORY argument is crucial to create an additional
# "smash" folder at destination!
install(DIRECTORY include/smash DESTINATION "include" PATTERN "*.h.in" EXCLUDE)
install(DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/include/smash" DESTINATION "include")
# Copy SMASH executable to installation bin folder
install(PROGRAMS $<TARGET_FILE_NAME:smash> DESTINATION bin)
