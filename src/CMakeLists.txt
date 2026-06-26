########################################################
#
#    Copyright (c) 2012-2026
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
add_executable(smash smash.cc)
add_library(smash_shared SHARED $<TARGET_OBJECTS:smash_objlib>)
add_library(smash_api INTERFACE)

# Alias for the shared library to unify build and install trees syntax. If e.g. SMASH is used via
# add_subdirectory in another project, the alias allows to link to SMASH always using SMASH::smash.
add_library(SMASH::smash ALIAS smash_shared)

# Both smash executable and shared library are linked privately to the object library because the
# object is only an implementation detail and we do not want to expose it to clients, but also to
# avoid that the dependency on the external libraries is propagated to users of the shared library,
# which should not be aware of them (dependency on external library is introduced in the Config
# file). What has to be propagated to clients belongs to smash_api interface target.
target_link_libraries(smash_objlib PUBLIC smash_project_build_requirements smash_api)
target_link_libraries(smash_api
                      INTERFACE smash_project_usage_requirements
                                einhard
                                yaml-cpp
                                cuhre
                                suave
                                divonne
                                vegas)
target_link_libraries(smash_shared PRIVATE smash_objlib PUBLIC smash_api)
target_link_libraries(smash PRIVATE smash_shared)

# Disable compiler-specific extensions for all targets. This cannot be done through INTERFACE
# libraries because properties like CXX_EXTENSIONS do not propagate through them.
set_target_properties(smash_objlib smash smash_shared PROPERTIES CXX_EXTENSIONS OFF)

# Call shared library smash (e.g. libsmash.so). Note that, since "the build-tree runtime path will
# not be used for binaries in the install tree" (CMake documentation), we request here to set the
# INSTALL_RPATH in the library to ensure that projects linking to SMASH finds at run time libraries
# like Pythia, which are likely to be in not standard places.
set_target_properties(smash_shared
                      PROPERTIES OUTPUT_NAME smash
                                 EXPORT_NAME smash
                                 INSTALL_RPATH_USE_LINK_PATH TRUE
                                 INSTALL_RPATH "$<IF:$<PLATFORM_ID:Darwin>,@loader_path,$ORIGIN>")
# For the executable, we set the RPATH to be the same as the shared library one, but with an
# additional "../lib" to look for the shared library in the installation lib folder. Note that
# INSTALL_RPATH_USE_LINK_PATH is set to TRUE for the executable as well, since it is needed to find
# at runtime the shared library dependencies of the executable (e.g. Pythia) as smash executable
# directly depends on e.g. Pythia and not only through the smash shared library.
set_target_properties(smash
                      PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE
                                 INSTALL_RPATH
                                 "$<IF:$<PLATFORM_ID:Darwin>,@loader_path/../lib,$ORIGIN/../lib>")

# Set up compilation include paths for given target. Note that we need to set both BUILD_INTERFACE
# and INSTALL_INTERFACE to be sure that the include directories are correctly set both when using
# the target from the build tree and from the install tree.
function(smash_set_target_include_directories target scope)
    target_include_directories(${target}
                               ${scope}
                               $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
                               $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
                               $<INSTALL_INTERFACE:include>)
endfunction()

# Set up compilation include paths for the API target
smash_set_target_include_directories(smash_api INTERFACE)

# Enable warnings for supported compilers only (we know the compilers support the flag, do not test)
if(CMAKE_CXX_COMPILER_ID MATCHES "^((Apple)?Clang|GNU)$")
    target_compile_options(smash_objlib
                           PUBLIC -W
                                  -Wall
                                  -Wextra
                                  -Wmissing-declarations
                                  -Wpointer-arith
                                  -Wshadow
                                  -Wuninitialized
                                  -Wfloat-conversion
                                  -Winit-self
                                  -Wundef
                                  -Wcast-align
                                  -Wformat=2
                                  -Wold-style-cast
                                  -Werror=switch)
    # GNU specific options
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        target_compile_options(smash_objlib
                               PUBLIC -Wlogical-op
                                      -Wnull-dereference
                                      -Wduplicated-cond
                                      -Wrestrict
                                      -Wduplicated-branches
                                      -Wsuggest-override)
        # Needed for workaround about C++17 compatibility (bug in some GNU versions)
        target_compile_definitions(smash_objlib PUBLIC GCC_COMPILER)
    else()
        # LLVM specific options
        target_compile_options(smash_objlib PUBLIC -Wno-missing-braces
                                                   -Winconsistent-missing-override)
        # Use of contracted mathematical expression like FMA (if the hardware supports them) are
        # controlled by a flag whose default unfortunately varies among compilers and versions
        # thereof. Since clang-14 changed the default from "off" to "on", we set here the value to
        # "on" also for older versions. GNU default is "fast", but this seems to be similar to
        # clang's "on" (which fuses in the same statement) than to clang's "fast" (which fuses
        # across statements). Here we simply try to have the same behaviour across the supported
        # compilers.
        if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 14)
            target_compile_options(smash_objlib PUBLIC -ffp-contract=on)
        endif()
    endif()
    # Enable warnings about deprecated C functions by default.
    option(DEPRECATE_C_FNS "Enable warnings about deprecated C functions in SMASH." OFF)
    if(DEPRECATE_C_FNS)
        target_compile_definitions(smash_objlib PUBLIC DONT_ASSERT_TEST)
        target_compile_options(smash_objlib PUBLIC -include smash/deprecate_c_functions.h)
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

# Note that we look for Pythia first (and especially before ROOT, HepMC3 and Rivet) to ensure it is
# linked before all the other and it is meant to avoid possible issues with e.g. ROOT whose include
# directories might contain a Pythia.h header that might not be the one of the Pythia linked to
# SMASH. See https://github.com/smash-transport/smash/issues/38#issuecomment-1033239643
set(Pythia_REQUIRED_VERSION "8.316") # For SMASHConfig.cmake file
find_package(Pythia ${Pythia_REQUIRED_VERSION} EXACT REQUIRED)
if(Pythia_FOUND)
    if(TARGET Pythia::Pythia)
        target_link_libraries(smash_api INTERFACE Pythia::Pythia)
        target_compile_definitions(smash_api INTERFACE "PYTHIA_XML_DIR=\"${Pythia_xmldoc_PATH}\""
                                                       "PYTHIA_FOUND")
    else()
        message(FATAL_ERROR "Pythia library found, but its target is somehow missing.\n"
                            "Ensure to be using the provided CMake module to find it.\n")
    endif()
endif()

set(GSL_MINIMUM_VERSION "2.0") # For SMASHConfig.cmake file
find_package(GSL ${GSL_MINIMUM_VERSION} REQUIRED)
if(GSL_FOUND)
    if(TARGET GSL::gsl AND TARGET GSL::gslcblas)
        target_link_libraries(smash_api INTERFACE GSL::gsl GSL::gslcblas)
    else()
        message(FATAL_ERROR " \n" " GSL library found, but its targets are somehow missing.\n"
                            " Ensure to be using the official CMake module to find it.\n")
    endif()
endif()

set(EIGEN3_MINIMUM_VERSION "3.0") # For SMASHConfig.cmake file
find_package(Eigen3 ${EIGEN3_MINIMUM_VERSION} REQUIRED)
if(EIGEN3_FOUND)
    if(TARGET Eigen3::Eigen)
        target_link_libraries(smash_api INTERFACE Eigen3::Eigen)
    else()
        message(FATAL_ERROR " \n" " Eigen3 library found, but its target is somehow missing.\n"
                            " Ensure to be using the provided CMake module to find it.\n")
    endif()
endif()

option(TRY_USE_ROOT "Turn this off to disable ROOT output support in SMASH." ON)
set(SMASH_USE_ROOT OFF) # For SMASHConfig.cmake file
set(ROOT_COMPONENTS Core RIO Tree) # Also for SMASHConfig.cmake file
if(TRY_USE_ROOT)
    find_package(ROOT CONFIG QUIET COMPONENTS ${ROOT_COMPONENTS})
    set(ROOT_MINIMUM_VERSION "6.20") # Also for SMASHConfig.cmake file
    if(ROOT_FOUND)
        message(STATUS "Candidate ROOT version ${ROOT_VERSION} found, "
                       "requested ${ROOT_MINIMUM_VERSION}")
        if(ROOT_VERSION VERSION_LESS ROOT_MINIMUM_VERSION)
            set(ROOT_FOUND FALSE)
        endif()
    endif()
    if(ROOT_FOUND)
        message(STATUS "Found ROOT: ${ROOT_INCLUDE_DIRS} (version ${ROOT_VERSION})")
        set(_root_targets ROOT::Core ROOT::RIO ROOT::Tree)
        # Prepare dedicated target for SMASH installation
        add_library(smash_root INTERFACE)
        add_library(SMASH::smash_root ALIAS smash_root)
        smash_set_target_include_directories(smash_root INTERFACE)
        target_compile_definitions(smash_root INTERFACE SMASH_USE_ROOT)
        target_link_libraries(smash_root INTERFACE ${_root_targets} SMASH::smash)
        set_target_properties(smash_root PROPERTIES EXPORT_NAME smash_root)
        set(SMASH_USE_ROOT ON)
        # Pull ROOT into smash_objlib to compile ROOT related translation units
        target_compile_definitions(smash_objlib PUBLIC SMASH_USE_ROOT)
        target_link_libraries(smash_objlib PUBLIC ${_root_targets})
        # Link ROOT into smash executable
        target_link_libraries(smash PRIVATE smash_root)
    else()
        message(STATUS "Suitable ROOT not found. Support disabled. "
                       "See README for more information.")
    endif()
endif()

option(TRY_USE_HEPMC "Turn this off to disable HepMC3 output support in SMASH." ON)
set(SMASH_USE_HepMC3 OFF) # For SMASHConfig.cmake file
if(TRY_USE_HEPMC OR TRY_USE_RIVET)
    find_package(HepMC3 QUIET CONFIG) # Min version not handled here but after to better inform user
    set(HepMC3_MINIMUM_VERSION "3.2.6") # Also for SMASHConfig.cmake file
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
            set(_hepmc_targets HepMC3::HepMC3)
            set(_hepmc_definitions SMASH_USE_HEPMC)
        endif()
        if(NOT ROOT_FOUND)
            message(STATUS "HepMC3_treeroot output disabled, missing ROOT support")
        elseif(NOT TARGET HepMC3::rootIO)
            message(STATUS "HepMC3_treeroot output disabled, missing HepMC3 ROOT I/O support")
        else()
            message(STATUS "HepMC3_treeroot output support enabled")
            list(APPEND _hepmc_targets HepMC3::rootIO)
            list(APPEND _hepmc_definitions SMASH_USE_HEPMC_ROOTIO)
        endif()
        # Prepare dedicated target for SMASH installation
        add_library(smash_hepmc3 INTERFACE)
        add_library(SMASH::smash_hepmc3 ALIAS smash_hepmc3)
        smash_set_target_include_directories(smash_hepmc3 INTERFACE)
        target_compile_definitions(smash_hepmc3 INTERFACE ${_hepmc_definitions})
        target_link_libraries(smash_hepmc3 INTERFACE ${_hepmc_targets} SMASH::smash)
        set_target_properties(smash_hepmc3 PROPERTIES EXPORT_NAME smash_hepmc3)
        set(SMASH_USE_HepMC3 ON)
        # Pull HepMC3 into smash_objlib to compile HepMC3 related translation units
        target_compile_definitions(smash_objlib PUBLIC ${_hepmc_definitions})
        target_link_libraries(smash_objlib PUBLIC ${_hepmc_targets})
        # Link HepMC3 into smash executable
        target_link_libraries(smash PRIVATE smash_hepmc3)
    else()
        message(STATUS "Suitable HepMC3 not found. Support disabled. "
                       "See README for more information.")
    endif()
endif()

option(TRY_USE_RIVET "Turn this off to disable Rivet output support in SMASH." ON)
set(SMASH_USE_Rivet OFF) # For SMASHConfig.cmake file
if(TRY_USE_RIVET AND HepMC3_FOUND)
    find_package(Rivet QUIET) # Min version not handled here but after to better inform user
    set(Rivet_MINIMUM_VERSION "3.1") # Also for SMASHConfig.cmake file
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
        # Prepare dedicated target for SMASH installation
        add_library(smash_rivet INTERFACE)
        add_library(SMASH::smash_rivet ALIAS smash_rivet)
        smash_set_target_include_directories(smash_rivet INTERFACE)
        target_compile_definitions(smash_rivet INTERFACE SMASH_USE_RIVET)
        target_link_libraries(smash_rivet INTERFACE Rivet::Rivet smash_hepmc3 SMASH::smash)
        set_target_properties(smash_rivet PROPERTIES EXPORT_NAME smash_rivet)
        set(SMASH_USE_Rivet ON)
        # Pull Rivet into smash_objlib to compile Rivet related translation units
        target_compile_definitions(smash_objlib PUBLIC SMASH_USE_RIVET)
        target_link_libraries(smash_objlib PUBLIC Rivet::Rivet)
        # Link Rivet into smash executable
        target_link_libraries(smash PRIVATE smash_rivet) # HepMC3 targets already linked in
    else()
        message(STATUS "Suitable Rivet not found. Support disabled. "
                       "See README for more information.")
    endif()
endif()

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
                       input_keys.cc
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

# Targets preparation for tests
if(BUILD_TESTING)
    # Add library for tests. Note that we link to smash_objlib publicly to propagate all SMASH
    # internal features to tests targets that will link against this library. Furthermore, we
    # exclude this library from all, as tests are and this library is to build tests only.
    add_library(smash_testlib EXCLUDE_FROM_ALL STATIC $<TARGET_OBJECTS:smash_objlib>)
    target_link_libraries(smash_testlib PUBLIC vir-test smash_objlib smash_api)
    set_target_properties(smash_testlib PROPERTIES OUTPUT_NAME smash-for-tests CXX_EXTENSIONS OFF)

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
    function(check_sanitizer_flag_support flag result_var)
        set(_saved_flags "${CMAKE_REQUIRED_FLAGS}")
        set(CMAKE_REQUIRED_FLAGS "${flag}") # add_compiler_flag only adds the flag for
        # compilation - it's also required with the link command, though
        add_compiler_flag("${flag}" CXX_RESULT _have_flag)
        set(CMAKE_REQUIRED_FLAGS "${_saved_flags}")
        set(${result_var} ${_have_flag} PARENT_SCOPE)
    endfunction()

    set(SANITIZER_FLAG "-fsanitize=address,undefined")
    check_sanitizer_flag_support("${SANITIZER_FLAG}" HAVE_SANITIZER)
    if(NOT HAVE_SANITIZER)
        set(SANITIZER_FLAG "-fsanitize=address")
        check_sanitizer_flag_support("${SANITIZER_FLAG}" HAVE_SANITIZER)
    endif()

    # Offer option to disable sanitizer usage, if available.
    if(HAVE_SANITIZER)
        option(USE_SANITIZER
               "Use available address and undefined behavior sanitizers for test code." ON)
    else()
        set(USE_SANITIZER FALSE)
    endif()

    # The static library for unit tests uses the same object files as smash_objlib but adds
    # test-specific definitions (BUILD_TESTS) and optional sanitizer instrumentation. Note that
    # smash_objlib itself is NOT modified, so main builds remain unaffected.
    target_compile_definitions(smash_testlib
                               PUBLIC "BUILD_TESTS"
                                      "TEST_CONFIG_PATH=std::filesystem::path(\"${PROJECT_SOURCE_DIR}\")"
    )

    if(USE_SANITIZER)
        target_compile_options(smash_testlib PUBLIC ${SANITIZER_FLAG})
        target_link_options(smash_testlib PUBLIC ${SANITIZER_FLAG})
        # Create a smash exectuable with sanitizers enabled to be used by tests. As for tests,
        # exclude this target from all target, since it is meant to be used by tests only.
        add_executable(smash_sanitizer EXCLUDE_FROM_ALL smash.cc)
        target_compile_definitions(smash_sanitizer PUBLIC BUILD_TESTS)
        target_compile_options(smash_sanitizer PUBLIC ${SANITIZER_FLAG})
        target_link_options(smash_sanitizer PUBLIC ${SANITIZER_FLAG})
        target_link_libraries(smash_sanitizer PRIVATE smash_testlib)
    endif()

    add_subdirectory(tests)
endif()

# Configure and install CMake configuration files for find_package.
configure_package_config_file(${CMAKE_SOURCE_DIR}/cmake/SMASHConfig.cmake.in
                              ${CMAKE_CURRENT_BINARY_DIR}/SMASHConfig.cmake
                              INSTALL_DESTINATION lib/cmake/SMASH)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/SMASHConfig.cmake DESTINATION lib/cmake/SMASH)

# Install targets and files. Note that we do not only install the shared library target, but also
# the API interface and the usage requirements interface, since this chain of dependencies will be
# needed by clients using the smash libraries in a CMake project.
install(TARGETS smash_shared smash_api smash_project_usage_requirements
        EXPORT SMASH-targets-core
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        RUNTIME DESTINATION bin)
# Note that the absent trailing slash in the DIRECTORY argument is crucial to create an additional
# "smash" folder at destination!
install(FILES ${generated_headers} DESTINATION "include/smash")
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/include/smash" DESTINATION include PATTERN "*.h.in"
                                                                                          EXCLUDE)
install(DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/include/smash" DESTINATION include)
# Install SMASH executable to installation bin folder
install(TARGETS smash RUNTIME DESTINATION bin)
# If optional dependencies are found, install them into dedicated targets files in order to split
# the installed SMASH into CMake components. This allows users in their project to e.g. use
# find_package(SMASH COMPONENTS Core ROOT) to get only the core and ROOT support without the need to
# link to HepMC3 and Rivet. The components names are defined in the SMASHConfig.cmake.in file.
foreach(_dependency IN ITEMS ROOT HepMC3 Rivet)
    string(TOLOWER ${_dependency} _dependency_lower)
    if(SMASH_USE_${_dependency})
        install(EXPORT SMASH-targets-${_dependency_lower}
                FILE SMASH-targets-${_dependency_lower}.cmake
                NAMESPACE SMASH::
                DESTINATION lib/cmake/SMASH)
        install(TARGETS smash_${_dependency_lower}
                EXPORT SMASH-targets-${_dependency_lower}
                INCLUDES
                DESTINATION include)
        if(_dependency STREQUAL "Rivet")
            install(FILES ${CMAKE_SOURCE_DIR}/cmake/FindRivet.cmake DESTINATION lib/cmake/SMASH)
        endif()
    endif()
endforeach()
