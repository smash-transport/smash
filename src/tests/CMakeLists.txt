########################################################
#
#    Copyright (c) 2014-2015
#      SMASH Team
#
#    BSD 3-clause license
#
#########################################################

# remove -Winline for test code
string(REPLACE "-Winline " "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")

file(COPY ../../input/config.yaml DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

string(REPLACE " -Wfloat-equal" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")

option(FORMAT_TEST_OUTPUT_FOR_VIM
   "If ON the unit test output of run_* targets will work with the vim quickfix parser and vim will directly take you to the line that failed."
   OFF)
set(UNITTEST_ARGS "")
if(FORMAT_TEST_OUTPUT_FOR_VIM)
   set(UNITTEST_ARGS "-v")
endif()


# tests:
# Helper macro that adds a run_<name> target
macro(smash_add_test name depends)
   add_test(NAME ${name} COMMAND ${ARGN})
   add_custom_target(run_${name}
      COMMAND ${ARGN} ${UNITTEST_ARGS}
      DEPENDS ${depends}
      COMMENT "Executing test ${name}"
      VERBATIM
      )
endmacro()

macro(smash_add_exe name)
  add_executable(${name} ${name}.cc)
  target_link_libraries(${name} smash_lib ${SMASH_LIBRARIES})
  if(USE_SANITIZER)
      set_target_properties(${name} PROPERTIES
         COMPILE_FLAGS "${SANITIZER_FLAG} -DSMASH_TEST_OUTPUT_PATH=\\\"${PROJECT_BINARY_DIR}/test_output/${name}\\\""
         LINK_FLAGS ${SANITIZER_FLAG}
         )
   else()
      set_target_properties(${name} PROPERTIES
         COMPILE_FLAGS "-DSMASH_TEST_OUTPUT_PATH=\\\"${PROJECT_BINARY_DIR}/test_output/${name}\\\""
         )
   endif()
endmacro()

macro(smash_add_unittest name)
   smash_add_exe(${name})
   smash_add_test(${name} ${name} ${name})
endmacro()

# helper macro to run smash with certain arguments, choosing an output directory according to the test name
macro(smash_add_runtest name depends)
   smash_add_test(${name} ${depends} ${ARGN} -o "${PROJECT_BINARY_DIR}/test_output/${name}" -f)
endmacro()

# test source code to be clean, requires cpplint in the path
find_program(CPPLINT cpplint.py)
if(NOT CPPLINT)
   find_program(CPPLINT cpplint)
endif()
if(CPPLINT)
   message(STATUS "cpplint.py found at ${CPPLINT}")
   set(CPPLINT_ARGS --filter=-build/c++11,-build/include_subdir,+build/include_alpha,-runtime/operator,-readability/streams,-readability/function,-runtime/references,-whitespace/operators --counting=detailed --repository=${PROJECT_SOURCE_DIR})
   file(GLOB CPPLINT_FILES
        ${PROJECT_SOURCE_DIR}/src/include/*.h ${PROJECT_SOURCE_DIR}/src/*.cc)
   add_custom_target(cpplint
                     COMMAND ${CPPLINT} ${CPPLINT_ARGS} ${CPPLINT_FILES}
                     WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                     COMMENT "Executing test cpplint")
   # make sure we always run cpplint as a test
   add_test(NAME cpplint
            COMMAND ${CPPLINT} ${CPPLINT_ARGS} ${CPPLINT_FILES}
            WORKINNG_DIRECTORY ${PROJECT_SOURCE_DIR})
   add_custom_target(cpplint_build
                     COMMAND ${CPPLINT} ${CPPLINT_ARGS} --filter=-build/c++11,-build/include,-readability,-runtime,-whitespace ${CPPLINT_FILES}
                     WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                     COMMENT "Executing test cpplint_build")
   add_custom_target(cpplint_readability
                     COMMAND ${CPPLINT} ${CPPLINT_ARGS} --filter=-build,-readability/streams,-runtime,-whitespace ${CPPLINT_FILES}
                     WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                     COMMENT "Executing test cpplint_readability")
   add_custom_target(cpplint_runtime
                     COMMAND ${CPPLINT} ${CPPLINT_ARGS} --filter=-build,-readability,-runtime/references,-whitespace ${CPPLINT_FILES}
                     WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                     COMMENT "Executing test cpplint_runtime")
   add_custom_target(cpplint_whitespace
                     COMMAND ${CPPLINT} ${CPPLINT_ARGS} --filter=-build,-readability,-runtime ${CPPLINT_FILES}
                     WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                     COMMENT "Executing test cpplint_whitespace")
   file(GLOB_RECURSE CPPLINT_FILES_FULL
        ${PROJECT_SOURCE_DIR}/src/*.h ${PROJECT_SOURCE_DIR}/src/*.cc)
   add_custom_target(cpplint_full
                     COMMAND ${CPPLINT} ${CPPLINT_ARGS} ${CPPLINT_FILES_FULL}
                     WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                     COMMENT "Executing test cpplint_full")
else()
   message(WARNING "cpplint.py was not found. The cpplint test will be disabled.\nDownload from: http://google-styleguide.googlecode.com/svn/trunk/cpplint/cpplint.py")
endif()

# static code analysis with cppcheck
find_program(CPPCHECK cppcheck)
if(CPPCHECK)
   message(STATUS "cppcheck found at ${CPPCHECK}")
   set(CPPCHECK_ARGS -q -U_LIBCPP_BEGIN_NAMESPACE_STD -U_LIBCPP_END_NAMESPACE_STD -I src/include/ -i src/tests/ src)
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
else()
   message(WARNING "cppcheck was not found. The cppcheck test will be disabled.")
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
   message(WARNING "clang-tidy was not found. The clang-tidy test will be disabled.")
endif()

# generate git statistics, requires gitstats in the path
find_program(GITSTATS gitstats)
if(GITSTATS)
   message(STATUS "gitstats found at ${GITSTATS}")
   add_custom_target(stats
                     COMMAND ${GITSTATS} ${PROJECT_SOURCE_DIR} ${PROJECT_BINARY_DIR}/gitstats
                     COMMENT "Generating git statistics")
else()
   message(WARNING "gitstats was not found. The gitstats support will be disabled.")
endif()

add_definitions("-DTEST_CONFIG_PATH=bf::path(\"${PROJECT_SOURCE_DIR}\")")
add_definitions("-DBUILD_TESTS")

# compile-only tests
smash_add_exe(angles_distribution)
smash_add_exe(angles_zero)
smash_add_exe(woods-saxon)

# unit tests for classes:
smash_add_unittest(action)
smash_add_unittest(actions)
smash_add_unittest(adaptiverejectionsampler)
smash_add_unittest(angles)
smash_add_unittest(average)
smash_add_unittest(binaryoutput)
smash_add_unittest(clebschgordan)
smash_add_unittest(clock)
smash_add_unittest(configuration)
smash_add_unittest(decayaction)
smash_add_unittest(decaymodes)
smash_add_unittest(decaytree)
smash_add_unittest(deformednucleus)
smash_add_unittest(density)
smash_add_unittest(dileptons)
smash_add_unittest(distributions)
smash_add_unittest(enable_float_traps)
smash_add_unittest(energymomentumtensor)
smash_add_unittest(experiment)
smash_add_unittest(filelock)
smash_add_unittest(fourvector)
smash_add_unittest(grandcan_thermalizer)
smash_add_unittest(grid)
smash_add_unittest(hadgas_eos)
smash_add_unittest(hadgas_eos2)
smash_add_unittest(initial_conditions)
smash_add_unittest(integrate)
smash_add_unittest(interpolation)
smash_add_unittest(isospin)
smash_add_unittest(kinematics)
smash_add_unittest(lattice)
smash_add_unittest(listmodus)
smash_add_unittest(lorentzboost)
smash_add_unittest(lowess)
smash_add_unittest(mass_sampling)
smash_add_unittest(nucleus)
smash_add_unittest(oscar2013output)
smash_add_unittest(oscar1999output)
smash_add_unittest(parametrizations)
smash_add_unittest(particledata)
smash_add_unittest(particles)
smash_add_unittest(particletype)
smash_add_unittest(pauliblocking)
smash_add_unittest(pdgcode)
smash_add_unittest(photons)
smash_add_unittest(potentials)
smash_add_unittest(processbranch)
smash_add_unittest(propagate)
smash_add_unittest(quantumnumbers)
smash_add_unittest(random)
smash_add_unittest(sanity)
smash_add_unittest(scatteraction)
smash_add_unittest(scatteractionsfinder)
smash_add_unittest(spectral_functions)
smash_add_unittest(stringfunctions)
smash_add_unittest(tabulation)
smash_add_unittest(threevector)
smash_add_unittest(two_unstable_products)
smash_add_unittest(vtkoutput)
smash_add_unittest(width)
smash_add_unittest(without_float_traps)
smash_add_unittest(yamltest)


# verify that the binary has a cli help
smash_add_runtest(smash_help smash smash -h)

# verify that the binary runs with the various modi:
smash_add_runtest(smash_run smash smash)
smash_add_runtest(box_run smash smash -m Box -e 100.0
                  -c "Modi: {Box: {Initial_Condition: thermal momenta}}"
                  -c "Modi: {Box: {Length: 10.0}}"
                  -c "Modi: {Box: {Temperature: 0.2}}"
                  -c "Modi: {Box: {Start_Time: 0.0}}"
                  -c "Modi: {Box: {Init_Multiplicities: {211: 100}}}"
                  -c "Modi: {Box: {Init_Multiplicities: {-211: 100}}}"
                  -c "Modi: {Box: {Init_Multiplicities: {111: 100}}}"
                  -c "Modi: {Box: {Init_Multiplicities: {2212: 50}}}"
                  -c "Modi: {Box: {Init_Multiplicities: {2112: 50}}}")
smash_add_runtest(collider_run smash smash -m Collider)
smash_add_runtest(collider_frozen_fermi_run smash smash -m Collider
                  -c "Modi: {Collider: {Fermi_Motion: frozen}}")
smash_add_runtest(sphere_run smash smash -m Sphere -e 100.0
                  -c "Modi: {Sphere: {Radius: 5.0}}"
                  -c "Modi: {Sphere: {Sphere_Temperature: 0.2}}"
                  -c "Modi: {Sphere: {Start_Time: 0.0}}"
                  -c "Modi: {Box: {Initial_Condition: thermal momenta}}"
                  -c "Modi: {Sphere: {Init_Multiplicities: {211: 100}}}"
                  -c "Modi: {Sphere: {Init_Multiplicities: {-211: 100}}}"
                  -c "Modi: {Sphere: {Init_Multiplicities: {111: 100}}}"
                  -c "Modi: {Sphere: {Init_Multiplicities: {2212: 50}}}"
                  -c "Modi: {Sphere: {Init_Multiplicities: {2112: 50}}}")
smash_add_runtest(dilepton_run smash smash -d ${CMAKE_CURRENT_SOURCE_DIR}/decaymodes_incl_dileptons.txt
                  -c "Output: {Dileptons: {Format: [Oscar2013]}}")

# verify that binary runs certain time (this corresponds to 20 steps)
smash_add_runtest(box_steps smash smash -e 2.0 -c "General: {Delta_Time: 0.1}"
                  -m Box -c "Modi: {Box: {Initial_Condition: thermal momenta}}"
                  -c "Modi: {Box: {Length: 10.0}}"
                  -c "Modi: {Box: {Temperature: 0.2}}"
                  -c "Modi: {Box: {Start_Time: 0.0}}"
                  -c "Modi: {Box: {Init_Multiplicities: {211: 100}}}"
                  -c "Modi: {Box: {Init_Multiplicities: {-211: 100}}}"
                  -c "Modi: {Box: {Init_Multiplicities: {111: 100}}}"
                  -c "Modi: {Box: {Init_Multiplicities: {2212: 50}}}"
                  -c "Modi: {Box: {Init_Multiplicities: {2112: 50}}}")
smash_add_runtest(collider_steps smash smash -e 2.0
                  -c "General: {Delta_Time: 0.1}" -m Collider)
smash_add_runtest(sphere_steps smash smash -e 2.0
                  -c "General: {Delta_Time: 0.1}" -m Sphere
                  -c "Modi: {Sphere: {Radius: 5.0}}"
                  -c "Modi: {Sphere: {Sphere_Temperature: 0.2}}"
                  -c "Modi: {Sphere: {Start_Time: 0.0}}"
                  -c "Modi: {Box: {Initial_Condition: thermal momenta}}"
                  -c "Modi: {Sphere: {Init_Multiplicities: {211: 100}}}"
                  -c "Modi: {Sphere: {Init_Multiplicities: {-211: 100}}}"
                  -c "Modi: {Sphere: {Init_Multiplicities: {111: 100}}}"
                  -c "Modi: {Sphere: {Init_Multiplicities: {2212: 50}}}"
                  -c "Modi: {Sphere: {Init_Multiplicities: {2112: 50}}}")

# verify that default run has no mem leaks
find_program(CTEST_MEMORYCHECK_COMMAND valgrind)
if(CTEST_MEMORYCHECK_COMMAND)
   message(STATUS "valgrind found at ${CTEST_MEMORYCHECK_COMMAND}")
   set(KNOWN_BUG "\n       Note: Known bug with valgrind 3.11 - unrecognised instruction.")
   set(CTEST_MEMORYCHECK_COMMAND_OPTIONS "-v" "--error-exitcode=1")
   add_custom_target(memcheck_collider
                     COMMAND ${CTEST_MEMORYCHECK_COMMAND} ${CTEST_MEMORYCHECK_COMMAND_OPTIONS}
                             ${PROJECT_BINARY_DIR}/smash
                             -c "General: {End_Time: 5.0, Delta_Time: 0.1}"
                             -o "${PROJECT_BINARY_DIR}/test_output/memcheck_collider"
                             -f
                     DEPENDS smash
                     COMMENT "Executing test memcheck_collider. ${KNOWN_BUG}"
                     VERBATIM)
   add_custom_target(memcheck_box
                     COMMAND ${CTEST_MEMORYCHECK_COMMAND} ${CTEST_MEMORYCHECK_COMMAND_OPTIONS}
                             ${PROJECT_BINARY_DIR}/smash
                             -c "General: {End_Time: 5.0, Delta_Time: 0.1}"
                             -i ${PROJECT_SOURCE_DIR}/input/box/config.yaml
                             -p ${PROJECT_SOURCE_DIR}/input/box/particles.txt
                             -d ${PROJECT_SOURCE_DIR}/input/box/decaymodes.txt
                             -o "${PROJECT_BINARY_DIR}/test_output/memcheck_box"
                             -f
                     DEPENDS smash
                     COMMENT "Executing test memcheck_box. ${KNOWN_BUG}"
                     VERBATIM)
   add_custom_target(memcheck_sphere
                     COMMAND ${CTEST_MEMORYCHECK_COMMAND} ${CTEST_MEMORYCHECK_COMMAND_OPTIONS}
                     ${PROJECT_BINARY_DIR}/smash
                     -c "General: {End_Time: 5.0, Delta_Time: 0.1}"
                     -i ${PROJECT_SOURCE_DIR}/input/sphere/config.yaml
                     -o "${PROJECT_BINARY_DIR}/test_output/memcheck_sphere" -f
                     DEPENDS smash
                     COMMENT "Executing test memcheck_sphere. ${KNOWN_BUG}"
                     VERBATIM)
else()
   message(STATUS "valgrind not found. Memcheck tests are disabled!")
endif()

if(USE_SANITIZER)
   smash_add_runtest(sanitizer_box smash_sanitizer smash_sanitizer -e 5.0
                    -m Box -c "Modi: {Box: {Initial_Condition: thermal momenta}}"
                    -c "Modi: {Box: {Length: 10.0}}"
                    -c "Modi: {Box: {Temperature: 0.2}}"
                    -c "Modi: {Box: {Start_Time: 0.0}}"
                    -c "Modi: {Box: {Init_Multiplicities: {211: 100}}}"
                    -c "Modi: {Box: {Init_Multiplicities: {-211: 100}}}"
                    -c "Modi: {Box: {Init_Multiplicities: {111: 100}}}"
                    -c "Modi: {Box: {Init_Multiplicities: {2212: 50}}}"
                    -c "Modi: {Box: {Init_Multiplicities: {2112: 50}}}"
                    -c "Collision_Term: {Strings: false}")
   smash_add_runtest(sanitizer_collider smash_sanitizer smash_sanitizer -e 5.0 -m Collider
                    -c "Collision_Term: {Strings: false}")
   smash_add_runtest(sanitizer_sphere   smash_sanitizer smash_sanitizer -e 5.0
                    -m Sphere -c "Modi: {Sphere: {Radius: 5.0}}"
                    -c "Modi: {Sphere: {Sphere_Temperature: 0.2}}"
                    -c "Modi: {Sphere: {Start_Time: 0.0}}"
                    -c "Modi: {Box: {Initial_Condition: thermal momenta}}"
                    -c "Modi: {Sphere: {Init_Multiplicities: {211: 100}}}"
                    -c "Modi: {Sphere: {Init_Multiplicities: {-211: 100}}}"
                    -c "Modi: {Sphere: {Init_Multiplicities: {111: 100}}}"
                    -c "Modi: {Sphere: {Init_Multiplicities: {2212: 50}}}"
                    -c "Modi: {Sphere: {Init_Multiplicities: {2112: 50}}}"
                    -c "Collision_Term: {Strings: false}")
endif()
