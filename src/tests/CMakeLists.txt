########################################################
#
#    Copyright (c) 2012-2026
#      SMASH Team
#
#    BSD 3-clause license
#
########################################################

# Copy needed input file to build directory
configure_file(../../input/config.yaml ${CMAKE_CURRENT_BINARY_DIR}/config.yaml COPYONLY)

# Add option for Vim users to format the unit test output in a way that works with the Vim quickfix
# parser. This allows them to directly jump to the line that failed in the test when running the
# run_* targets. This serves few developers, but the maintenance cost is very low.
option(FORMAT_TEST_OUTPUT_FOR_VIM
       "If ON the unit test output of run_* targets will work with the vim quickfix parser and vim will directly take you to the line that failed."
       OFF)

# Load functionality to add different types of tests and setup the environment for functional tests
include(TestsUtilities)

# -------------------------------- Unit tests  --------------------------------
smash_add_unittest(action)
smash_add_unittest(actions)
smash_add_unittest(alphaclusterednucleus)
smash_add_unittest(angles)
smash_add_unittest(angles_distribution)
smash_add_unittest(angles_zero)
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
smash_add_unittest(woods-saxon)
smash_add_unittest(yamltest)

# ----------------------------- Functional tests  -----------------------------
option(ENABLE_FUNCTIONAL_TESTS "Turn this on to enable SMASH functional tests." OFF)
if(ENABLE_FUNCTIONAL_TESTS)
    do_setup_for_functional_tests(setup_was_successful)
    if(setup_was_successful)
        smash_add_functionaltest(ic_for_hybrid)
    endif()
else()
    message(STATUS "Functional tests disabled.")
endif()

# -------------------------------- Run tests  ---------------------------------
smash_add_runtest(smash_help smash smash -h)
# Verify that the binary runs with the various modi. Here we use the shipped default configurations
# file, while particles and decay modes files are taken from the input directory.
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
smash_add_runtest(potentials_run smash smash -i ${PROJECT_SOURCE_DIR}/input/potentials/config.yaml)
smash_add_runtest(deformednucleus_run smash smash
                  -i ${PROJECT_SOURCE_DIR}/input/deformed_nucleus/config.yaml)
# Further verify that smash binary runs for a certain time (this corresponds to 20 steps)
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
# Verify that also other behavior of smash executable runs without crashing
smash_add_runtest(particle_table_dump smash smash -x)
smash_add_runtest(reactions_dump smash smash -l)
smash_add_runtest(crosssection_dump smash smash -s 2112,2212)
smash_add_runtest(resonance_spectral_function_dump smash smash -r 113)
# If sanitizer is enabled, repeat few runs using it
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
# If ROOT, HepMC3 and Rivet are found, also make a run using them (through output files)
if(ROOT_FOUND AND HepMC3_FOUND AND Rivet_FOUND)
    smash_add_runtest(third_party_output_run smash smash
                      -i ${PROJECT_SOURCE_DIR}/input/config.yaml
                      -c "Output: {Particles: {Format: [Root, HepMC_treeroot, HepMC_asciiv3]}}"
                      -c "Output: {Collisions: {Format: [Root, HepMC_treeroot, HepMC_asciiv3]}}"
                      -c "Output: {Rivet: {Format: [YODA]}}"
                      -c "Output: {Rivet: {Analyses: [MC_FSPARTICLES]}}")
endif()
# Finally mimic an external project using SMASH as a library, verifying that few executables linking
# to smash library can be built and run without crashing.
add_test(NAME usage_of_SMASH_as_library
         COMMAND bash run_examples_for_testing.bash "${PROJECT_SOURCE_DIR}"
                 "${Pythia_CONFIG_EXECUTABLE}"
         WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/examples/using_SMASH_as_library)
set_tests_properties(usage_of_SMASH_as_library PROPERTIES ENVIRONMENT
                                                          SMASH_BUILD_DIR=${PROJECT_BINARY_DIR})

# =================================================================================================

# Helper function that calls a command with --version option and save the output in a variable
function(get_version_string command output_variable)
    execute_process(COMMAND ${command} --version OUTPUT_VARIABLE result
                    ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)
    set(${output_variable} "${result}" PARENT_SCOPE)
endfunction()

# Setup tests to check that default run has no memory leaks, if valgrind is available
find_program(VALGRIND_EXECUTABLE valgrind)
if(VALGRIND_EXECUTABLE)
    message(STATUS "Valgrind found at ${VALGRIND_EXECUTABLE}")
    get_version_string(${VALGRIND_EXECUTABLE} VALGRIND_VERSION_STRING)
    if(VALGRIND_VERSION_STRING MATCHES "^valgrind-3.11")
        message(ATTENTION "There is a known bug with valgrind 3.11"
                          " (unrecognised instruction) - disabling memchecks.")
    else()
        set(CTEST_MEMORYCHECK_COMMAND ${VALGRIND_EXECUTABLE})
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

# Variable defining the regular expression for a generic version number
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
