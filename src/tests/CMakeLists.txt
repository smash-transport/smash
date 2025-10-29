########################################################
#
#    Copyright (c) 2012-2024
#      SMASH Team
#
#    BSD 3-clause license
#
########################################################

# Add 3rd-party virtest library
include_directories(SYSTEM "${CMAKE_SOURCE_DIR}/3rdparty/virtest")

# Copy needed input file to build directory
file(COPY ../../input/config.yaml DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

option(FORMAT_TEST_OUTPUT_FOR_VIM
       "If ON the unit test output of run_* targets will work with the vim quickfix parser and vim will directly take you to the line that failed."
       OFF)
set(UNITTEST_ARGS "")
if(FORMAT_TEST_OUTPUT_FOR_VIM)
    set(UNITTEST_ARGS "-v")
endif()

# tests: Helper macro that adds a run_<name> target
macro(smash_add_test name depends)
    add_test(NAME ${name} COMMAND ${ARGN})
    if(APPLE)
        # Known false positive container overflow error of AddressSanitizer on Mac.
        # https://github.com/google/sanitizers/wiki/AddressSanitizerContainerOverflow
        set_property(TEST ${name} PROPERTY ENVIRONMENT "ASAN_OPTIONS=detect_container_overflow=0")
    endif()
    add_custom_target(run_${name}
                      COMMAND ${ARGN} ${UNITTEST_ARGS}
                      DEPENDS ${depends}
                      COMMENT "Executing test ${name}"
                      VERBATIM)
endmacro()

macro(smash_add_exe name)
    add_executable(${name} ${name}.cc)
    target_link_libraries(${name} smash_static ${SMASH_LIBRARIES})
    if(USE_SANITIZER)
        set_target_properties(${name}
                              PROPERTIES COMPILE_FLAGS
                                         "${SANITIZER_FLAG} -DSMASH_TEST_OUTPUT_PATH=\\\"${PROJECT_BINARY_DIR}/test_output/${name}\\\""
                                         LINK_FLAGS ${SANITIZER_FLAG})
    else()
        set_target_properties(${name}
                              PROPERTIES COMPILE_FLAGS
                                         "-DSMASH_TEST_OUTPUT_PATH=\\\"${PROJECT_BINARY_DIR}/test_output/${name}\\\""
        )
    endif()
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        # Disable variable tracking in debug mode (never used so far and gcc gives annoying notes
        # when failing to track them)
        target_compile_options(${name} PRIVATE "-fno-var-tracking-assignments")
        # The GNU -Wduplicated-branches flag makes compilation time of the lattice.cc test file
        # increase by order(s) of magnitude and it cannot basically be used
        if(${name} STREQUAL "lattice")
            target_compile_options(${name} PRIVATE "-Wno-duplicated-branches")
        endif()
    endif()
    if(CMAKE_CXX_COMPILER_ID MATCHES "^((Apple)?Clang|GNU)$")
        # It turned out that the -ffp-contract option of GNU and LLVM compilers (whose default is
        # different in the two compilers and among versions) is affecting reproducibility of results
        # when e.g. optimizing for different architectures (e.g. -march=native VS -march=x86-64).
        # Therefore, we decided to switch this feature off for tests to ensure that tests behave in
        # the same way on all machines.
        target_compile_options(${name} PRIVATE "-ffp-contract=off")
    endif()
endmacro()

macro(smash_add_unittest name)
    smash_add_exe(${name})
    smash_add_test(${name} ${name} ${name})
endmacro()

macro(smash_add_functionaltest name)
    add_test(NAME functional_${name}
             COMMAND ${Python3_EXECUTABLE} "${PROJECT_SOURCE_DIR}/src/tests/functional/${name}.py"
                     --source "${PROJECT_SOURCE_DIR}" --binary "${PROJECT_BINARY_DIR}")
    set_tests_properties(functional_${name} PROPERTIES FIXTURES_REQUIRED fixture_compile_smash)
endmacro()

# helper macro to run smash with certain arguments, choosing an output directory according to the
# test name
macro(smash_add_runtest name depends)
    smash_add_test(${name}
                   ${depends}
                   ${ARGN}
                   -o
                   "${PROJECT_BINARY_DIR}/test_output/${name}"
                   -f)
endmacro()

# helper macro that calls a command with --version option and save the output in a variable
macro(get_version_string command output_variable)
    execute_process(COMMAND ${command} --version OUTPUT_VARIABLE ${output_variable}
                    ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)
endmacro()

# variable defining the regular expression for a generic version number
set(GENERIC_VERSION_REGEX "[0-9]+([.][0-9]+([.][0-9]+([.][0-9]+)?)?)?")

# test source code to be clean, requires cpplint in the path
set(CPPLINT_VERSION_REQUIRED "1.6.0")
find_program(CPPLINT cpplint NAMES cpplint cpplint.py)
if(CPPLINT)
    message(STATUS "cpplint found at ${CPPLINT}")
    get_version_string(${CPPLINT} CPPLINT_VERSION_STRING)
    string(REGEX MATCH "cpplint (${GENERIC_VERSION_REGEX})" CPPLINT_VERSION_STRING
                 "${CPPLINT_VERSION_STRING}")
    set(CPPLINT_VERSION_FOUND "${CMAKE_MATCH_1}")
    if(NOT CPPLINT_VERSION_FOUND MATCHES "^${GENERIC_VERSION_REGEX}$")
        message(ATTENTION "Unable to extract cpplint version - disabling cpplint checks.")
    elseif(NOT CPPLINT_VERSION_FOUND VERSION_EQUAL "${CPPLINT_VERSION_REQUIRED}")
        message(ATTENTION "cpplint version (${CPPLINT_VERSION_FOUND}) different "
                          "from required one (${CPPLINT_VERSION_REQUIRED})"
                          " - disabling cpplint checks.")
    else()
        set(CPPLINT_ARGS
            --filter=-build/c++11,-build/include_subdir,+build/include_alpha,-runtime/operator,-readability/streams,-readability/function,-runtime/references,-whitespace/operators
            --counting=detailed
            --repository=${PROJECT_SOURCE_DIR})
        file(GLOB CPPLINT_FILES ${PROJECT_SOURCE_DIR}/src/include/smash/*.h
             ${PROJECT_SOURCE_DIR}/src/*.cc)
        add_custom_target(cpplint
                          COMMAND ${CPPLINT} ${CPPLINT_ARGS} ${CPPLINT_FILES}
                          WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                          COMMENT "Executing test cpplint")
        # make sure we always run cpplint as a test
        add_test(NAME cpplint COMMAND ${CPPLINT} ${CPPLINT_ARGS} ${CPPLINT_FILES}
                 WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
        add_custom_target(cpplint_build
                          COMMAND ${CPPLINT} ${CPPLINT_ARGS}
                                  --filter=-build/c++11,-build/include,-readability,-runtime,-whitespace
                                  ${CPPLINT_FILES}
                          WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                          COMMENT "Executing test cpplint_build")
        add_custom_target(cpplint_readability
                          COMMAND ${CPPLINT} ${CPPLINT_ARGS}
                                  --filter=-build,-readability/streams,-runtime,-whitespace
                                  ${CPPLINT_FILES}
                          WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                          COMMENT "Executing test cpplint_readability")
        add_custom_target(cpplint_runtime
                          COMMAND ${CPPLINT} ${CPPLINT_ARGS}
                                  --filter=-build,-readability,-runtime/references,-whitespace
                                  ${CPPLINT_FILES}
                          WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                          COMMENT "Executing test cpplint_runtime")
        add_custom_target(cpplint_whitespace
                          COMMAND ${CPPLINT} ${CPPLINT_ARGS} --filter=-build,-readability,-runtime
                                  ${CPPLINT_FILES}
                          WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                          COMMENT "Executing test cpplint_whitespace")
        file(GLOB_RECURSE CPPLINT_FILES_FULL ${PROJECT_SOURCE_DIR}/src/include/smash/*.h
             ${PROJECT_SOURCE_DIR}/src/*.cc)
        add_custom_target(cpplint_full
                          COMMAND ${CPPLINT} ${CPPLINT_ARGS} ${CPPLINT_FILES_FULL}
                          WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                          COMMENT "Executing test cpplint_full")
    endif()
else()
    message(ATTENTION "cpplint not found (version ${CPPLINT_VERSION_REQUIRED}"
                      " needed). The cpplint tests will be disabled.")
endif()

# static code analysis with cppcheck
set(CPPCHECK_VERSION_REQUIRED "2.8")
find_program(CPPCHECK cppcheck)
if(CPPCHECK)
    message(STATUS "cppcheck found at ${CPPCHECK}")
    get_version_string(${CPPCHECK} CPPCHECK_VERSION_STRING)
    string(REGEX MATCH "^Cppcheck (${GENERIC_VERSION_REGEX})$" CPPCHECK_VERSION_STRING
                 "${CPPCHECK_VERSION_STRING}")
    set(CPPCHECK_VERSION_FOUND "${CMAKE_MATCH_1}")
    if(NOT CPPCHECK_VERSION_FOUND MATCHES "^${GENERIC_VERSION_REGEX}$")
        message(ATTENTION "Unable to extract cppcheck version - disabling cppcheck checks.")
    elseif(NOT CPPCHECK_VERSION_FOUND VERSION_EQUAL "${CPPCHECK_VERSION_REQUIRED}")
        message(ATTENTION "cppcheck version (${CPPCHECK_VERSION_FOUND}) different "
                          "from required one (${CPPCHECK_VERSION_REQUIRED})"
                          " - disabling cppcheck tests.")
    else()
        set(CPPCHECK_ARGS
            -q
            -U_LIBCPP_BEGIN_NAMESPACE_STD
            -U_LIBCPP_END_NAMESPACE_STD
            -I
            src/include/
            -i
            src/tests/
            src)
        add_custom_target(cppcheck
                          COMMAND ${CPPCHECK} "--enable=all" ${CPPCHECK_ARGS}
                          WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                          COMMENT "Executing test cppcheck")
        add_custom_target(cppcheck_error
                          COMMAND ${CPPCHECK} ${CPPCHECK_ARGS}
                          WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                          COMMENT "Executing test cppcheck_error")
        add_custom_target(cppcheck_warning
                          COMMAND ${CPPCHECK} "--enable=warning" ${CPPCHECK_ARGS}
                          WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                          COMMENT "Executing test cppcheck_warning")
        add_custom_target(cppcheck_performance
                          COMMAND ${CPPCHECK} "--enable=performance" ${CPPCHECK_ARGS}
                          WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                          COMMENT "Executing test cppcheck_performance")
        add_custom_target(cppcheck_portability
                          COMMAND ${CPPCHECK} "--enable=portability" ${CPPCHECK_ARGS}
                          WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                          COMMENT "Executing test cppcheck_portability")
        add_custom_target(cppcheck_information
                          COMMAND ${CPPCHECK} "--enable=information" ${CPPCHECK_ARGS}
                          WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                          COMMENT "Executing test cppcheck_information")
        add_custom_target(cppcheck_style
                          COMMAND ${CPPCHECK} "--enable=style" ${CPPCHECK_ARGS}
                          WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                          COMMENT "Executing test cppcheck_style")
        add_custom_target(cppcheck_unused
                          COMMAND ${CPPCHECK} "--enable=unusedFunction" ${CPPCHECK_ARGS}
                          WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                          COMMENT "Executing test cppcheck_unused")
    endif()
else()
    message(ATTENTION "cppcheck not found (version ${CPPCHECK_VERSION_REQUIRED}"
                      " needed). The cppcheck tests will be disabled.")
endif()

# one more linter tool: clang-tidy
find_program(CLANGTIDY clang-tidy NAMES clang-tidy clang-tidy-3.6 clang-tidy-3.5)
if(CLANGTIDY)
    message(STATUS "clang-tidy found at ${CLANGTIDY}")
    file(GLOB CLANGTIDY_FILES ${PROJECT_SOURCE_DIR}/src/*.cc)
    add_custom_target(clang-tidy
                      COMMAND ${CLANGTIDY} ${CLANGTIDY_FILES} "-p=${PROJECT_BINARY_DIR}"
                      WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/src/
                      COMMENT "Executing test clang-tidy")
else()
    message(ATTENTION "clang-tidy not found. The clang-tidy test will be disabled.")
endif()

# generate git statistics, requires gitstats in the path
find_program(GITSTATS gitstats)
if(GITSTATS)
    message(STATUS "gitstats found at ${GITSTATS}")
    add_custom_target(stats COMMAND ${GITSTATS} ${PROJECT_SOURCE_DIR} ${PROJECT_BINARY_DIR}/gitstats
                      COMMENT "Generating git statistics")
else()
    message(ATTENTION "gitstats not found. The gitstats target will not be created.")
endif()

add_definitions("-DTEST_CONFIG_PATH=std::filesystem::path(\"${PROJECT_SOURCE_DIR}\")")
add_definitions("-DBUILD_TESTS")

# compile-only tests
smash_add_exe(angles_distribution)
smash_add_exe(angles_zero)
smash_add_exe(woods-saxon)

# unit tests for classes:
smash_add_unittest(action)
smash_add_unittest(actions)
smash_add_unittest(alphaclusterednucleus)
smash_add_unittest(angles)
smash_add_unittest(average)
smash_add_unittest(binaryoutput)
smash_add_unittest(clebschgordan)
smash_add_unittest(clebschgordan_lookup)
smash_add_unittest(clock)
smash_add_unittest(configuration)
smash_add_unittest(decayaction)
smash_add_unittest(decaymodes)
smash_add_unittest(decaytree)
smash_add_unittest(deformednucleus)
smash_add_unittest(density)
smash_add_unittest(dileptons)
smash_add_unittest(distributions)
smash_add_unittest(dynamic_fluidization)
smash_add_unittest(enable_float_traps)
smash_add_unittest(energymomentumtensor)
smash_add_unittest(experiment)
smash_add_unittest(filelock)
smash_add_unittest(formfactors)
smash_add_unittest(fourvector)
smash_add_unittest(icoutput)
smash_add_unittest(grandcan_thermalizer)
smash_add_unittest(grid)
smash_add_unittest(hadgas_eos)
smash_add_unittest(hadgas_eos2)
smash_add_unittest(hypersurfacecrossing)
smash_add_unittest(initial_conditions)
smash_add_unittest(input_keys)
smash_add_unittest(integrate)
smash_add_unittest(interpolation)
smash_add_unittest(interpolation2D)
smash_add_unittest(isospin)
smash_add_unittest(key)
smash_add_unittest(kinematics)
smash_add_unittest(lattice)
smash_add_unittest(listmodus)
smash_add_unittest(lorentzboost)
smash_add_unittest(lowess)
smash_add_unittest(mass_sampling)
smash_add_unittest(nucleus)
smash_add_unittest(numeric_cast)
smash_add_unittest(oscar2013output)
smash_add_unittest(oscar1999output)
smash_add_unittest(outputformatter)
smash_add_unittest(parametrizations)
smash_add_unittest(particledata)
smash_add_unittest(particles)
smash_add_unittest(particletype)
smash_add_unittest(pauliblocking)
smash_add_unittest(pdgcode)
smash_add_unittest(photons)
smash_add_unittest(potentials)
smash_add_unittest(processbranch)
smash_add_unittest(stringprocess)
smash_add_unittest(propagate)
smash_add_unittest(quantumnumbers)
smash_add_unittest(quantumsampling)
smash_add_unittest(random)
smash_add_unittest(rootsolver)
smash_add_unittest(sanity)
smash_add_unittest(scatteraction)
smash_add_unittest(scatteractionmulti)
smash_add_unittest(scatteractionsfinder)
smash_add_unittest(sha256)
smash_add_unittest(spheremodus)
smash_add_unittest(spectral_functions)
smash_add_unittest(stringfunctions)
smash_add_unittest(tabulation)
smash_add_unittest(threevector)
smash_add_unittest(traits)
smash_add_unittest(two_unstable_products)
smash_add_unittest(vtkoutput)
smash_add_unittest(width)
smash_add_unittest(without_float_traps)
smash_add_unittest(yamltest)

# verify that the binary has a cli help
smash_add_runtest(smash_help smash smash -h)

# Verify that the binary runs with the various modi; reading the shipped default configs, particles
# and decaymodes files from the input directory:
smash_add_runtest(smash_run smash smash)
smash_add_runtest(box_run smash smash
                  -i ${PROJECT_SOURCE_DIR}/input/box/config.yaml
                  -p ${PROJECT_SOURCE_DIR}/input/box/particles.txt
                  -d ${PROJECT_SOURCE_DIR}/input/box/decaymodes.txt)
smash_add_runtest(stochastic_box_run smash smash
                  -i ${PROJECT_SOURCE_DIR}/input/stochastic_box/config.yaml
                  -p ${PROJECT_SOURCE_DIR}/input/stochastic_box/particles_only_pi0.txt
                  -d ${PROJECT_SOURCE_DIR}/input/stochastic_box/decaymodes_all_off.txt)
smash_add_runtest(multi_particle_box_run smash smash
                  -i ${PROJECT_SOURCE_DIR}/input/multi_particle_box/config.yaml
                  -p ${PROJECT_SOURCE_DIR}/input/multi_particle_box/particles.txt
                  -d ${PROJECT_SOURCE_DIR}/input/multi_particle_box/decaymodes.txt)
smash_add_runtest(collider_run smash smash -m Collider)
smash_add_runtest(collider_frozen_fermi_run smash smash -m Collider
                  -c "Modi: {Collider: {Fermi_Motion: frozen}}")
smash_add_runtest(sphere_run smash smash -i ${PROJECT_SOURCE_DIR}/input/sphere/config.yaml)
smash_add_runtest(list_run smash smash -i ${PROJECT_SOURCE_DIR}/input/list/config.yaml
                  -c "Modi: {List: {File_Directory: ${PROJECT_SOURCE_DIR}/input/list}}")
smash_add_runtest(dilepton_run smash smash -i ${PROJECT_SOURCE_DIR}/input/dileptons/config.yaml
                  -d ${PROJECT_SOURCE_DIR}/input/dileptons/decaymodes.txt)
smash_add_runtest(light_nuclei_run smash smash
                  -i ${PROJECT_SOURCE_DIR}/input/config.yaml
                  -p ${PROJECT_SOURCE_DIR}/input/light_nuclei/particles.txt
                  -d ${PROJECT_SOURCE_DIR}/input/light_nuclei/decaymodes.txt)
smash_add_runtest(custom_nucleus_run smash smash
                  -i ${PROJECT_SOURCE_DIR}/input/custom_nucleus/config.yaml
                  -c "General: {Nevents: 2}"
                  -c "Modi: {Collider: {Projectile: {Custom: {File_Directory: ${PROJECT_SOURCE_DIR}/input/custom_nucleus}}}}"
                  -c "Modi: {Collider: {Target: {Custom: {File_Directory: ${PROJECT_SOURCE_DIR}/input/custom_nucleus}}}}"
)

# Setup of functional tests. Since these are written in Python, here the interpreter is found, a
# virtual environment is setup, the interpreter of the virtual environment is looked for and the
# needed Python requirements are installed. If any step fails, functional tests are disabled.
option(ENABLE_FUNCTIONAL_TESTS "Turn this on to enable SMASH functional tests." OFF)
if(ENABLE_FUNCTIONAL_TESTS)
    message(STATUS "Looking for Python3")
    find_package(Python3 3.3 QUIET COMPONENTS Interpreter Development)
    if(Python3_FOUND)
        # Set up virtual environment and skip functional tests if anything fails
        set(PYTHON_VENV_PATH "${PROJECT_BINARY_DIR}/python_venv")
        execute_process(COMMAND ${Python3_EXECUTABLE} -m venv "${PYTHON_VENV_PATH}"
                        RESULT_VARIABLE venv_result)
        if(venv_result AND NOT venv_result EQUAL 0)
            message(ATTENTION "Setting up Python3 venv failed. Functional tests disabled.")
            set(ENABLE_FUNCTIONAL_TESTS OFF)
        endif()
        # Look again for Python, this time to locate that of the venv (taken from CMake discourse)
        set(SYSTEM_PYTHON_EXECUTABLE ${Python3_EXECUTABLE})
        set(ENV{VIRTUAL_ENV} "${PYTHON_VENV_PATH}") # Update the environment with VIRTUAL_ENV
                                                    # variable (mimic the activate script)
        set(Python3_FIND_VIRTUALENV FIRST) # Change the context of the search
        unset(Python3_EXECUTABLE) # Needed as it is also an input variable for find_package
        find_package(Python3 QUIET COMPONENTS Interpreter Development)
        if(SYSTEM_PYTHON_EXECUTABLE STREQUAL Python3_EXECUTABLE)
            message(FATAL_ERROR " \n" " Python3 venv executable is the same as the system one.\n"
                                " This should not happen and needs to be investigated.\n")
        endif()
        message(STATUS "Python3 venv for functional tests setup at '${PYTHON_VENV_PATH}'.")
        message(STATUS "Found Python3: ${Python3_EXECUTABLE}")
    else()
        message(ATTENTION "Python3 not found. Functional tests disabled.")
        set(ENABLE_FUNCTIONAL_TESTS OFF)
    endif()
endif()
if(ENABLE_FUNCTIONAL_TESTS)
    message(STATUS "Installing Python3 requirements for functional tests into venv")
    execute_process(COMMAND ${Python3_EXECUTABLE} -m pip install --upgrade pip -r
                            "${PROJECT_SOURCE_DIR}/src/tests/functional/requirements.txt"
                    RESULT_VARIABLE pip_install_result
                    OUTPUT_FILE "${PYTHON_VENV_PATH}/pip_install_requirements.log")
    if(pip_install_result AND NOT pip_install_result EQUAL 0)
        message(ATTENTION "Installing Python3 requirements failed. Functional tests disabled.")
        set(ENABLE_FUNCTIONAL_TESTS OFF)
    else()
        message(STATUS "Installing Python3 requirements for functional tests into venv - done")
    endif()
endif()
# Finally create functional tests. Note that it is not necessary to activate or deactivate the
# virtual environment because we run functional tests explicitly using the Python_EXECUTABLE which
# is now set to that of the virtual environment.
if(ENABLE_FUNCTIONAL_TESTS)
    # Since all functional tests need smash executable to be compiled before running, we want ctest
    # to first compile SMASH and then run the functional tests. In order to do so we use tests
    # fixtures setting up a "test" to compile SMASH and using this as setup fixture. See CMake doc
    # https://cmake.org/cmake/help/latest/prop_test/FIXTURES_REQUIRED.html for further information.
    # Please note that we do not want to always compile SMASH before each test and hence we wrap the
    # SMASH compilation into a small bash test to skip it if the smash exectuable exists.
    set(BASH_CODE_FOR_TEST
        "[[ ! -f '${CMAKE_BINARY_DIR}/smash' ]] && \
            ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --config ${CMAKE_BUILD_TYPE} \
                             --target smash  || exit 0")
    add_test(NAME compile_smash_executable COMMAND bash -c ${BASH_CODE_FOR_TEST})
    set_tests_properties(compile_smash_executable PROPERTIES FIXTURES_SETUP fixture_compile_smash)
    smash_add_functionaltest(ic_for_hybrid)
else()
    message(STATUS "Functional tests disabled.")
endif()

if(ROOT_FOUND AND HepMC3_FOUND AND Rivet_FOUND)
    smash_add_runtest(third_party_output_run smash smash
                      -i ${PROJECT_SOURCE_DIR}/input/config.yaml
                      -c "Output: {Particles: {Format: [Root, HepMC_treeroot, HepMC_asciiv3]}}"
                      -c "Output: {Collisions: {Format: [Root, HepMC_treeroot, HepMC_asciiv3]}}"
                      -c "Output: {Rivet: {Format: [YODA]}}"
                      -c "Output: {Rivet: {Analyses: [MC_FSPARTICLES]}}")
endif()

# Test the shipped config files for potentials and deformed nuclei by verifying the binary runs with
# them.
smash_add_runtest(potentials_run smash smash -i ${PROJECT_SOURCE_DIR}/input/potentials/config.yaml)
smash_add_runtest(deformednucleus_run smash smash
                  -i ${PROJECT_SOURCE_DIR}/input/deformed_nucleus/config.yaml)

# verify that binary runs certain time (this corresponds to 20 steps)
smash_add_runtest(box_steps smash smash
                  -e 2.0
                  -c "General: {Delta_Time: 0.1}"
                  -m Box
                  -c "Modi: {Box: {Initial_Condition: thermal momenta}}"
                  -c "Modi: {Box: {Length: 10.0}}"
                  -c "Modi: {Box: {Temperature: 0.2}}"
                  -c "Modi: {Box: {Start_Time: 0.0}}"
                  -c "Modi: {Box: {Init_Multiplicities: {211: 100}}}"
                  -c "Modi: {Box: {Init_Multiplicities: {-211: 100}}}"
                  -c "Modi: {Box: {Init_Multiplicities: {111: 100}}}"
                  -c "Modi: {Box: {Init_Multiplicities: {2212: 50}}}"
                  -c "Modi: {Box: {Init_Multiplicities: {2112: 50}}}")
smash_add_runtest(collider_steps smash smash
                  -e 2.0
                  -c "General: {Delta_Time: 0.1}"
                  -m Collider)
smash_add_runtest(sphere_steps smash smash
                  -e 2.0
                  -c "General: {Delta_Time: 0.1}"
                  -m Sphere
                  -c "Modi: {Sphere: {Radius: 5.0}}"
                  -c "Modi: {Sphere: {Temperature: 0.2}}"
                  -c "Modi: {Sphere: {Start_Time: 0.0}}"
                  -c "Modi: {Box: {Initial_Condition: thermal momenta}}"
                  -c "Modi: {Sphere: {Init_Multiplicities: {211: 100}}}"
                  -c "Modi: {Sphere: {Init_Multiplicities: {-211: 100}}}"
                  -c "Modi: {Sphere: {Init_Multiplicities: {111: 100}}}"
                  -c "Modi: {Sphere: {Init_Multiplicities: {2212: 50}}}"
                  -c "Modi: {Sphere: {Init_Multiplicities: {2112: 50}}}")
smash_add_runtest(particle_table_dump smash smash -x)
smash_add_runtest(reactions_dump smash smash -l)
smash_add_runtest(crosssection_dump smash smash -s 2112,2212)
smash_add_runtest(resonance_spectral_function_dump smash smash -r 113)

# Setup tests to check that default run has no memory leaks, if valgrind is available
find_program(CTEST_MEMORYCHECK_COMMAND valgrind)
if(CTEST_MEMORYCHECK_COMMAND)
    message(STATUS "Valgrind found at ${CTEST_MEMORYCHECK_COMMAND}")
    get_version_string(${CTEST_MEMORYCHECK_COMMAND} VALGRIND_VERSION_STRING)
    if(VALGRIND_VERSION_STRING MATCHES "^valgrind-3.11")
        message(ATTENTION "There is a known bug with valgrind 3.11"
                          " (unrecognised instruction) - disabling memchecks.")
    else()
        set(CTEST_MEMORYCHECK_COMMAND_OPTIONS "-v" "--error-exitcode=1")
        add_custom_target(memcheck_collider
                          COMMAND ${CTEST_MEMORYCHECK_COMMAND} ${CTEST_MEMORYCHECK_COMMAND_OPTIONS}
                                  ${PROJECT_BINARY_DIR}/smash -c
                                  "General: {End_Time: 5.0, Delta_Time: 0.1}" -o
                                  "${PROJECT_BINARY_DIR}/test_output/memcheck_collider" -f
                          DEPENDS smash
                          COMMENT "Executing test memcheck_collider."
                          VERBATIM)
        add_custom_target(memcheck_box
                          COMMAND ${CTEST_MEMORYCHECK_COMMAND} ${CTEST_MEMORYCHECK_COMMAND_OPTIONS}
                                  ${PROJECT_BINARY_DIR}/smash -c
                                  "General: {End_Time: 5.0, Delta_Time: 0.1}" -i
                                  ${PROJECT_SOURCE_DIR}/input/box/config.yaml -p
                                  ${PROJECT_SOURCE_DIR}/input/box/particles.txt -d
                                  ${PROJECT_SOURCE_DIR}/input/box/decaymodes.txt -o
                                  "${PROJECT_BINARY_DIR}/test_output/memcheck_box" -f
                          DEPENDS smash
                          COMMENT "Executing test memcheck_box."
                          VERBATIM)
        add_custom_target(memcheck_sphere
                          COMMAND ${CTEST_MEMORYCHECK_COMMAND} ${CTEST_MEMORYCHECK_COMMAND_OPTIONS}
                                  ${PROJECT_BINARY_DIR}/smash -c
                                  "General: {End_Time: 5.0, Delta_Time: 0.1}" -i
                                  ${PROJECT_SOURCE_DIR}/input/sphere/config.yaml -o
                                  "${PROJECT_BINARY_DIR}/test_output/memcheck_sphere" -f
                          DEPENDS smash
                          COMMENT "Executing test memcheck_sphere."
                          VERBATIM)
    endif()
else()
    message(STATUS "Valgrind not found. Memcheck tests are disabled!")
endif()

add_test(NAME usage_of_SMASH_as_library
         COMMAND bash run_examples_for_testing.bash "${PROJECT_SOURCE_DIR}"
                 "${Pythia_CONFIG_EXECUTABLE}"
         WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/examples/using_SMASH_as_library)

set_tests_properties(usage_of_SMASH_as_library PROPERTIES ENVIRONMENT
                                                          SMASH_BUILD_DIR=${PROJECT_BINARY_DIR})

if(USE_SANITIZER)
    smash_add_runtest(sanitizer_box smash_sanitizer smash_sanitizer
                      -e 5.0
                      -m Box
                      -c "Modi: {Box: {Initial_Condition: thermal momenta}}"
                      -c "Modi: {Box: {Length: 10.0}}"
                      -c "Modi: {Box: {Temperature: 0.2}}"
                      -c "Modi: {Box: {Start_Time: 0.0}}"
                      -c "Modi: {Box: {Init_Multiplicities: {211: 100}}}"
                      -c "Modi: {Box: {Init_Multiplicities: {-211: 100}}}"
                      -c "Modi: {Box: {Init_Multiplicities: {111: 100}}}"
                      -c "Modi: {Box: {Init_Multiplicities: {2212: 50}}}"
                      -c "Modi: {Box: {Init_Multiplicities: {2112: 50}}}"
                      -c "Collision_Term: {Strings: false}")
    smash_add_runtest(sanitizer_collider smash_sanitizer smash_sanitizer
                      -e 5.0
                      -m Collider
                      -c "Collision_Term: {Strings: false}")
    smash_add_runtest(sanitizer_sphere smash_sanitizer smash_sanitizer
                      -e 5.0
                      -m Sphere
                      -c "Modi: {Sphere: {Radius: 5.0}}"
                      -c "Modi: {Sphere: {Temperature: 0.2}}"
                      -c "Modi: {Sphere: {Start_Time: 0.0}}"
                      -c "Modi: {Box: {Initial_Condition: thermal momenta}}"
                      -c "Modi: {Sphere: {Init_Multiplicities: {211: 100}}}"
                      -c "Modi: {Sphere: {Init_Multiplicities: {-211: 100}}}"
                      -c "Modi: {Sphere: {Init_Multiplicities: {111: 100}}}"
                      -c "Modi: {Sphere: {Init_Multiplicities: {2212: 50}}}"
                      -c "Modi: {Sphere: {Init_Multiplicities: {2112: 50}}}"
                      -c "Collision_Term: {Strings: false}")
endif()
