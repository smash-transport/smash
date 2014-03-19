# remove -Winline for test code
string(REPLACE "-Winline " "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")

# test sources
add_executable(fourvector fourvector.cc)
add_executable(lorentzboost lorentzboost.cc ../fourvector.cc ../random.cc)
add_executable(angles_distribution angles_distribution.cc ../random.cc)
add_executable(angles_zero angles_zero.cc ../random.cc)
add_executable(clebschgordan clebschgordan.cc ../clebschgordan.cc)

file(COPY ../config_general.yaml DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
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

macro(smash_add_unittest name)
   add_executable(${name} ${name}.cc)
   set_target_properties(${name} PROPERTIES
      COMPILE_FLAGS ${ASAN_FLAG}
      LINK_FLAGS ${ASAN_FLAG}
      )
   target_link_libraries(${name} smash_lib ${SMASH_LIBRARIES})
   smash_add_test(${name} ${name} ${name})
endmacro()

# test source code to be clean, requires cpplint in the path
find_program(CPPLINT cpplint.py)
if(CPPLINT)
   file(GLOB_RECURSE CPPLINT_ARGS ${PROJECT_SOURCE_DIR}/src/*.h 
      ${PROJECT_SOURCE_DIR}/src/*.cc)
   set(CPPLINT_ARGS --filter=-build/include_order,-readability/streams ${CPPLINT_ARGS})
   message(STATUS "cpplint.py found at ${CPPLINT}")
   add_test(NAME cclint COMMAND ${CPPLINT} ${CPPLINT_ARGS}
      WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
      )
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

smash_add_unittest(yamltest)
smash_add_unittest(configuration)
smash_add_unittest(experiment)
smash_add_unittest(particles)
smash_add_unittest(decaymodes)

# some angle operations tests
smash_add_unittest(angles)
smash_add_unittest(nucleus)
smash_add_unittest(woods-saxon)
smash_add_unittest(random)

# some fourvector operations tests
smash_add_test(fourvector fourvector fourvector)

# Lorentz boost tests
smash_add_test(lorentzboost lorentzboost lorentzboost)

# clebsch-gordan coefficient tests
smash_add_test(clebschgordan clebschgordan clebschgordan)

set(SMASH_OPTIONS -o "${PROJECT_BINARY_DIR}/test_output" -f)

# verify that the binary has a cli help
smash_add_test(mash_help smash smash -h)

# verify that the binary runs
smash_add_test(mash_run smash smash ${SMASH_OPTIONS})

# verify that the binary can run a box calculation
smash_add_test(box_run smash smash -m Box ${SMASH_OPTIONS})

# verify that the binary can do Collider
smash_add_test(collider_run smash smash -m Collider ${SMASH_OPTIONS})

# verify that the binary can do sphere
#smash_add_test(sphere_run smash smash -m Sphere ${SMASH_OPTIONS})

# verify that binary runs certain number of steps
smash_add_test(box_steps smash smash -s 20 -m Box ${SMASH_OPTIONS})
smash_add_test(collider_steps smash smash -s 20 -m Collider ${SMASH_OPTIONS})
#smash_add_test(sphere_steps smash smash -s 20 -m Sphere ${SMASH_OPTIONS})

# verify that default run has no mem leaks
find_program(CTEST_MEMORYCHECK_COMMAND valgrind)
if(CTEST_MEMORYCHECK_COMMAND)
   message(STATUS "valgrind found at ${CTEST_MEMORYCHECK_COMMAND}")
   set(CTEST_MEMORYCHECK_COMMAND_OPTIONS "-v")
   smash_add_test(memcheck_box smash ${CTEST_MEMORYCHECK_COMMAND}
      ${CTEST_MEMORYCHECK_COMMAND_OPTIONS}
      ${PROJECT_BINARY_DIR}/smash -c "General: {STEPS: 10, MODUS: Box}" ${SMASH_OPTIONS})
   smash_add_test(memcheck_collider smash ${CTEST_MEMORYCHECK_COMMAND}
      ${CTEST_MEMORYCHECK_COMMAND_OPTIONS}
      ${PROJECT_BINARY_DIR}/smash -s 50 -m Collider ${SMASH_OPTIONS})
#  smash_add_test(memcheck_sphere smash ${CTEST_MEMORYCHECK_COMMAND}
#     ${CTEST_MEMORYCHECK_COMMAND_OPTIONS}
#     ${PROJECT_BINARY_DIR}/smash -s 50 -m Sphere ${SMASH_OPTIONS})
else()
   message(STATUS "valgrind not found. Memcheck tests are disabled!")
endif()

if(HAVE_ADDRESS_SANITIZER)
   smash_add_test(asan_box      smash_asan smash_asan -c "General: {STEPS: 50,MODUS: Box}" ${SMASH_OPTIONS})
   smash_add_test(asan_collider smash_asan smash_asan -s 500 -m Collider ${SMASH_OPTIONS})
   #smash_add_test(asan_sphere   smash_asan smash_asan -s 500 -m Sphere ${SMASH_OPTIONS})
endif()
