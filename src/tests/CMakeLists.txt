# remove -Winline for test code
string(REPLACE "-Winline " "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")

file(COPY ../config.yaml DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ../particles.txt DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ../decaymodes.txt DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

string(REPLACE " -Wfloat-equal" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")

# tests:
# Helper macro that adds a run_<name> target
macro(smash_add_test name depends)
   add_test(NAME ${name} COMMAND ${ARGN})
   add_custom_target(run_${name}
      COMMAND ${ARGN}
      DEPENDS ${depends}
      COMMENT "Executing test ${name}"
      VERBATIM
      )
endmacro()

macro(smash_add_exe name)
  add_executable(${name} ${name}.cc)
  target_link_libraries(${name} smash_lib ${SMASH_LIBRARIES})
   if(HAVE_SANITIZER)
      set_target_properties(${name} PROPERTIES
         COMPILE_FLAGS "${SANITIZER_FLAG} -DSMASH_TEST_OUTPUT_PATH=\\\"${PROJECT_BINARY_DIR}/${name}_output/\\\""
         LINK_FLAGS ${SANITIZER_FLAG}
         )
   else()
      set_target_properties(${name} PROPERTIES
         COMPILE_FLAGS "-DSMASH_TEST_OUTPUT_PATH=\\\"${PROJECT_BINARY_DIR}/${name}_output/\\\""
         )
   endif()
endmacro()

macro(smash_add_unittest name)
   smash_add_exe(${name})
   smash_add_test(${name} ${name} ${name})
endmacro()

# test source code to be clean, requires cpplint in the path
find_program(CPPLINT cpplint.py)
if(CPPLINT)
   file(GLOB_RECURSE CPPLINT_ARGS ${PROJECT_SOURCE_DIR}/src/*.h
      ${PROJECT_SOURCE_DIR}/src/*.cc)
   set(CPPLINT_ARGS --filter=-build/include_order,-readability/streams,-build/include ${CPPLINT_ARGS})
   message(STATUS "cpplint.py found at ${CPPLINT}")
   add_custom_target(run_cclint
      COMMAND ${CPPLINT} ${CPPLINT_ARGS}
      WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
      COMMENT "Executing test cclint"
      )
else()
   message(WARNING "cpplint.py was not found. The cclint test will be disabled.\nDownload from: http://google-styleguide.googlecode.com/svn/trunk/cpplint/cpplint.py")
endif()

add_definitions("-DTEST_CONFIG_PATH=\"${PROJECT_SOURCE_DIR}/src\"")
add_definitions("-DBUILD_TESTS")

# compile-only tests
smash_add_exe(angles_distribution)
smash_add_exe(angles_zero)
smash_add_exe(woods-saxon)

# unit tests for classes:
smash_add_unittest(angles)
smash_add_unittest(binaryoutput)
smash_add_unittest(clebschgordan)
smash_add_unittest(clock)
smash_add_unittest(configuration)
smash_add_unittest(decaymodes)
smash_add_unittest(deformednucleus)
smash_add_unittest(distributions)
smash_add_unittest(experiment)
smash_add_unittest(fourvector)
smash_add_unittest(initial_conditions)
smash_add_unittest(lorentzboost)
smash_add_unittest(nucleus)
smash_add_unittest(oscar2013output)
smash_add_unittest(oscaroutput)
smash_add_unittest(particledata)
smash_add_unittest(particles)
smash_add_unittest(particletype)
smash_add_unittest(pdgcode)
smash_add_unittest(processbranch)
smash_add_unittest(propagate)
smash_add_unittest(quantumnumbers)
smash_add_unittest(random)
smash_add_unittest(sanity)
smash_add_unittest(threevector)
smash_add_unittest(vtkoutput)
smash_add_unittest(yamltest)

# for now -o has to be "data", because OSCAR output still expects to output to
# $PWD/data. In the future we want to use "${PROJECT_BINARY_DIR}/test_output"
set(SMASH_OPTIONS -o "data" -f)

# verify that the binary has a cli help
smash_add_test(mash_help smash smash -h)

# verify that the binary runs with the various modi:
smash_add_test(mash_run smash smash ${SMASH_OPTIONS})
smash_add_test(box_run smash smash -m Box ${SMASH_OPTIONS})
smash_add_test(collider_run smash smash -m Collider ${SMASH_OPTIONS})
smash_add_test(nucleus_run smash smash -m Nucleus ${SMASH_OPTIONS})
smash_add_test(sphere_run smash smash -m Sphere ${SMASH_OPTIONS})

# verify that binary runs certain time (this corresponds to 20 steps)
smash_add_test(box_steps smash smash -e 2.0 -c "General: {DELTA_TIME: 0.1}" -m Box ${SMASH_OPTIONS})
smash_add_test(collider_steps smash smash -e 2.0 -c "General: {DELTA_TIME: 0.1}" -m Collider ${SMASH_OPTIONS})
smash_add_test(nucleus_steps smash smash -e 2.0 -c "General: {DELTA_TIME: 0.1}" -m Nucleus ${SMASH_OPTIONS})
smash_add_test(sphere_steps smash smash -e 2.0 -c "General: {DELTA_TIME: 0.1}" -m Sphere ${SMASH_OPTIONS})

# verify that default run has no mem leaks
find_program(CTEST_MEMORYCHECK_COMMAND valgrind)
if(CTEST_MEMORYCHECK_COMMAND)
   message(STATUS "valgrind found at ${CTEST_MEMORYCHECK_COMMAND}")
   set(CTEST_MEMORYCHECK_COMMAND_OPTIONS "-v")
   smash_add_test(memcheck_box smash ${CTEST_MEMORYCHECK_COMMAND}
      ${CTEST_MEMORYCHECK_COMMAND_OPTIONS}
      ${PROJECT_BINARY_DIR}/smash -c "General: {END_TIME: 1.0, MODUS: Box}" ${SMASH_OPTIONS})
   smash_add_test(memcheck_collider smash ${CTEST_MEMORYCHECK_COMMAND}
      ${CTEST_MEMORYCHECK_COMMAND_OPTIONS}
      ${PROJECT_BINARY_DIR}/smash -c "General: {END_TIME: 5.0, DELTA_TIME: 0.1}" -m Collider ${SMASH_OPTIONS})
   smash_add_test(memcheck_nucleus smash ${CTEST_MEMORYCHECK_COMMAND}
      ${CTEST_MEMORYCHECK_COMMAND_OPTIONS}
      ${PROJECT_BINARY_DIR}/smash -c "General: {END_TIME: 5.0, DELTA_TIME: 0.1}" -m Nucleus ${SMASH_OPTIONS})
   smash_add_test(memcheck_sphere smash ${CTEST_MEMORYCHECK_COMMAND}
      ${CTEST_MEMORYCHECK_COMMAND_OPTIONS}
      ${PROJECT_BINARY_DIR}/smash -c "General: {END_TIME: 5.0, DELTA_TIME: 0.1}" -m Sphere ${SMASH_OPTIONS})
else()
   message(STATUS "valgrind not found. Memcheck tests are disabled!")
endif()

if(HAVE_SANITIZER)
   smash_add_test(sanitizer_box      smash_sanitizer smash_sanitizer -e 5.0 -m Box ${SMASH_OPTIONS})
   smash_add_test(sanitizer_collider smash_sanitizer smash_sanitizer -e 5.0 -m Collider ${SMASH_OPTIONS})
   smash_add_test(sanitizer_nucleus  smash_sanitizer smash_sanitizer -e 5.0 -m Nucleus ${SMASH_OPTIONS})
   smash_add_test(sanitizer_sphere   smash_sanitizer smash_sanitizer -e 5.0 -m Sphere ${SMASH_OPTIONS})
endif()
