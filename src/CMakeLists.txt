########################################################
#
#    Copyright (c) 2012-2024
#      SMASH Team
#
#    BSD 3-clause license
#
#########################################################

include(CheckLibraryExists)
include(AddCompilerFlag)
include(FindPkgConfig)

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
    option(DEPRECATE_C_FNS "Turn this on to enable warning about deprecated C funtions in SMASH."
           OFF)
    if(DEPRECATE_C_FNS)
        add_definitions(-include smash/deprecate_c_functions.h -DDONT_ASSERT_TEST)
    endif()
endif()

find_package(Git)
if(GIT_FOUND)
    # message("git found: ${GIT_EXECUTABLE} in version ${GIT_VERSION_STRING}")
    get_filename_component(GIT_PARENT_DIR ${CMAKE_CURRENT_SOURCE_DIR} PATH)
    set(GIT_DIR "${GIT_PARENT_DIR}/.git")
    if(EXISTS "${GIT_DIR}")
        execute_process(COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD
                        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} OUTPUT_VARIABLE GIT_BRANCH)
        string(STRIP ${GIT_BRANCH} GIT_BRANCH)
    else()
        set(GIT_BRANCH "NOTFOUND")
    endif()
    # message("git branch: ${GIT_BRANCH}")
endif(GIT_FOUND)

string(TIMESTAMP BUILD_DATE)

find_package(GSL 2.0 REQUIRED)
find_package(Eigen3 3.0 REQUIRED)

option(TRY_USE_ROOT "Turn this off to disable ROOT output support in SMASH." ON)
if(TRY_USE_ROOT)
    find_package(ROOT QUIET) # Min version not handled here but after to better inform user
    set(ROOT_MINIMUM_VERSION "5.34")
    if(ROOT_FOUND)
        message(STATUS "Candidate ROOT version ${ROOT_VERSION} found, "
                       "requested ${ROOT_MINIMUM_VERSION}")
        if(${ROOT_VERSION} VERSION_LESS ${ROOT_MINIMUM_VERSION})
            set(ROOT_FOUND FALSE)
        endif()
    endif()
    if(ROOT_FOUND)
        message(STATUS "Found valid ROOT ${ROOT_VERSION} (include at ${ROOT_INCLUDE_DIR}).")
        include_directories(SYSTEM "${ROOT_INCLUDE_DIR}")
        set(SMASH_LIBRARIES ${ROOT_LIBRARIES})
        add_definitions(-DSMASH_USE_ROOT)
    else()
        message(STATUS "Suitable ROOT not found. Support disabled. "
                       "See README for more information.")
    endif()
endif()

option(TRY_USE_HEPMC "Turn this off to disable HepMC output support in SMASH." ON)
if(TRY_USE_HEPMC OR TRY_USE_RIVET)
    find_package(HepMC3 QUIET) # Min version not handled here but after to better inform user
    set(HepMC3_MINIMUM_VERSION "3.2.3")
    if(HepMC3_FOUND)
        message(STATUS "Candidate HepMC version ${HepMC3_VERSION} found, "
                       "requested ${HepMC3_MINIMUM_VERSION}")
        if(${HepMC3_VERSION} VERSION_LESS ${HepMC3_MINIMUM_VERSION})
            set(HepMC3_FOUND FALSE)
        endif()
    endif()
    if(HepMC3_FOUND)
        message(STATUS "Found valid HepMC ${HepMC3_VERSION} (include at ${HEPMC3_INCLUDE_DIR}).")
        include_directories(SYSTEM "${HEPMC3_INCLUDE_DIR}")
        set(SMASH_LIBRARIES ${SMASH_LIBRARIES} ${HEPMC3_LIBRARIES})
        add_definitions(-DSMASH_USE_HEPMC)
        if((NOT HEPMC3_ROOTIO_LIB) OR (NOT ROOT_FOUND))
            message(STATUS "HepMC3_treeroot output disabled, missing ROOT and/or "
                           "HepMC3 ROOT I/O support")
        else()
            message(STATUS "HepMC3_treeroot output support enabled")
            add_definitions(-DSMASH_USE_HEPMC_ROOTIO)
        endif()
    else()
        message(STATUS "Suitable HepMC not found. Support disabled. "
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
        if(${RIVET_VERSION} VERSION_LESS "${Rivet_MINIMUM_VERSION}")
            set(Rivet_FOUND FALSE)
        elseif(${RIVET_VERSION} VERSION_GREATER_EQUAL "3.2")
            message(ATTENTION "Rivet version larger than 3.1.x found. "
                              "The SMASH interface might miss some functionality.")
        endif()
    endif()
    if(Rivet_FOUND)
        message(STATUS "Found valid Rivet ${RIVET_VERSION} (include at ${RIVET_INCLUDE_DIR}).")
        include_directories(SYSTEM "${RIVET_INCLUDE_DIR}")
        set(SMASH_LIBRARIES ${SMASH_LIBRARIES} ${RIVET_LIBRARIES})
        add_definitions(-DSMASH_USE_RIVET)
        set(Rivet_FOUND TRUE)
    else()
        message(STATUS "Suitable Rivet not found. Support disabled. "
                       "See README for more information.")
    endif()
endif()

# find Pythia
find_package(Pythia 8.315 EXACT REQUIRED)
if(Pythia_FOUND)
    # the BEFORE at the line below is due to
    # https://github.com/smash-transport/smash/issues/38#issuecomment-1033239643
    include_directories(BEFORE SYSTEM "${Pythia_INCLUDE_DIRS}")
    set(SMASH_LIBRARIES ${SMASH_LIBRARIES} ${Pythia_LIBRARIES} ${CMAKE_DL_LIBS})
    add_definitions("-DPYTHIA_XML_DIR=\"${Pythia_xmldoc_PATH}\"" -DPYTHIA_FOUND)
endif()

# set up include paths
include_directories(include)
include_directories("${CMAKE_CURRENT_BINARY_DIR}/include")
include_directories(SYSTEM ${GSL_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIR})

# set default libraries for linking
set(SMASH_LIBRARIES
    ${SMASH_LIBRARIES}
    ${GSL_LIBRARY}
    ${GSL_CBLAS_LIBRARY}
    einhard
    yaml-cpp
    cuhre
    suave
    divonne
    vegas # Cuba multidimensional integration
)

# ~~~
# For GNU and Clang compilers, link to C++ filesystem library if compiler version requires it
#  -> See "Notes" at https://en.cppreference.com/w/cpp/filesystem
#
# NOTE: As of Xcode 11 and MacOS 10.15, <filesystem> is shipped in the main
#       libc++.dylib and there's no need to link any additional library.
# ~~~
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9.1)
        set(SMASH_LIBRARIES ${SMASH_LIBRARIES} stdc++fs)
    endif()
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9.0)
        set(SMASH_LIBRARIES ${SMASH_LIBRARIES} c++fs)
    endif()
endif()

# list the source files
set(smash_src
    action.cc
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
    set(smash_src ${smash_src} rootoutput.cc)
endif()

if(TRY_USE_HEPMC AND HepMC3_FOUND)
    set(smash_src ${smash_src} hepmcoutput.cc hepmcinterface.cc)
    if(TRY_USE_RIVET AND Rivet_FOUND)
        set(smash_src ${smash_src} rivetoutput.cc)
    endif()
endif()

option(ENABLE_NANOBENCHMARKING "Turn this on to enable code to perform nanobenchmarking in SMASH."
       OFF)
if(ENABLE_NANOBENCHMARKING)
    if(CMAKE_HOST_SYSTEM_PROCESSOR MATCHES "^arm")
        message(WARNING "Nanobenchmarking is not supported on ARM processors and it was disabled.")
    else()
        set(smash_src ${smash_src} tsc.cc)
        add_definitions(-DUSE_NANOBENCHMARKING_CODE)
    endif()
endif()

# this is the "object library" target: compiles the sources only once see
# https://stackoverflow.com/a/29824424 NOTE: shared libraries need PIC and this is already set for
# all targets through CMAKE_POSITION_INDEPENDENT_CODE
add_library(objlib OBJECT ${smash_src})

add_executable(smash smash.cc $<TARGET_OBJECTS:objlib>)

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

    # define BUILD_TESTS for smash_static in order to enable small code changes for the unit tests.
    # It is an unfortunate side effect that the sanitzer build now also compiles to slightly
    # different code. To fix this we'd have to compile all source files a third time... not sure we
    # prefer that.
    if(USE_SANITIZER)
        set_target_properties(smash_static
                              PROPERTIES COMPILE_FLAGS "${SANITIZER_FLAG} -DBUILD_TESTS"
                                         LINK_FLAGS ${SANITIZER_FLAG})

        add_executable(smash_sanitizer smash.cc)
        set_target_properties(smash_sanitizer
                              PROPERTIES COMPILE_FLAGS "${SANITIZER_FLAG} -DBUILD_TESTS"
                                         LINK_FLAGS ${SANITIZER_FLAG})
        target_link_libraries(smash_sanitizer smash_static ${SMASH_LIBRARIES})
    else()
        set_target_properties(smash_static PROPERTIES COMPILE_FLAGS "-DBUILD_TESTS")
    endif()

    add_subdirectory(tests)
endif()

# Since at the top level the dependency of the install target on all target was disabled, it is
# needed to now build the smash library and exectuable target here before installing them
install(CODE "execute_process(COMMAND ${CMAKE_COMMAND}
                                      --build ${PROJECT_BINARY_DIR}
                                      --target smash_shared smash)")
install(TARGETS smash_shared LIBRARY DESTINATION "lib/${SMASH_INSTALLATION_SUBFOLDER}")
install(FILES ${generated_headers} DESTINATION "include/${SMASH_INSTALLATION_SUBFOLDER}/smash")
# Note that the absent trailing slash in the DIRECTORY argument is crucial to create an additional
# "smash" folder at destination!
install(DIRECTORY include/smash DESTINATION "include/${SMASH_INSTALLATION_SUBFOLDER}"
        PATTERN "*.h.in" EXCLUDE)
install(DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/include/smash"
        DESTINATION "include/${SMASH_INSTALLATION_SUBFOLDER}")
# Copy SMASH executable to installation bin folder adding a suffix to allow different versions to
# coexist. A symbolic link called 'smash' is also created, pointing to the last installed version.
# Due to this mechanism, we do not use "install(TARGETS ...)", as it does not offer a renaming
# mechanism.
install(PROGRAMS $<TARGET_FILE_NAME:smash> DESTINATION bin RENAME "${SMASH_INSTALLATION_SUBFOLDER}")
install(CODE "execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
                                      ${SMASH_INSTALLATION_SUBFOLDER}
                                      smash
                              WORKING_DIRECTORY ${CMAKE_INSTALL_PREFIX}/bin)")
