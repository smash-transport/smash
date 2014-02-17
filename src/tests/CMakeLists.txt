# test sources
add_executable(fourvector fourvector.cc)
add_executable(lorentzboost lorentzboost.cc ../fourvector.cc)
add_executable(angles_distribution angles_distribution.cc)
add_executable(particles particles.cc ../particles.cc ../fourvector.cc ../distributions.cc ../outputroutines.cc ../time.cc)
add_executable(clebschgordan clebschgordan.cc ../clebschgordan.cc)

file(COPY ../config_general.txt DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ../particles.txt DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ../decaymodes.txt DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

target_link_libraries(particles ${SMASH_LIBRARIES})

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
   message(STATUS "cpplint.py found at ${CPPLINT}")
   #file(GLOB source "../*.cc" "../include/*.h" "*.cc")
   add_test(NAME cclint COMMAND ${CPPLINT} */*.cc */*/*.h
      WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
      )
   add_custom_target(run_cclint
      COMMAND ${CPPLINT} */*.cc */*/*.cc */*/*.h
      WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
      COMMENT "Executing test cclint"
      )
else()
   message(WARNING "cpplint.py was not found. The cclint test will be disabled.\nDownload from: http://google-styleguide.googlecode.com/svn/trunk/cpplint/cpplint.py")
endif()

smash_add_unittest(experiment)

# some angle operations tests
smash_add_unittest(angles)

# some fourvector operations tests
smash_add_test(fourvector fourvector fourvector)

# Lorentz boost tests
smash_add_test(lorentzboost lorentzboost lorentzboost)


# some particles operations tests
smash_add_test(particles particles particles)

# clebsch-gordan coefficient tests
smash_add_test(clebschgordan clebschgordan clebschgordan)

# verify that the binary has a cli help
smash_add_test(mash_help smash smash -h)

# verify that the binary runs
smash_add_test(mash_run smash smash)

# verify that the binary can run a box calculation
smash_add_test(box_run smash smash -m Box)

# verify that the binary can do Collider
smash_add_test(collider_run smash smash -m Collider)

# verify that the binary can do sphere
#smash_add_test(sphere_run smash smash -m Sphere)

# verify that binary runs certain number of steps
smash_add_test(box_steps smash smash -s 20 -m Box)
smash_add_test(collider_steps smash smash -s 20 -m Collider)
#smash_add_test(sphere_steps smash smash -s 20 -m Sphere)

# verify that default run has no mem leaks
find_program(CTEST_MEMORYCHECK_COMMAND valgrind)
if(CTEST_MEMORYCHECK_COMMAND)
   message(STATUS "valgrind found at ${CTEST_MEMORYCHECK_COMMAND}")
   set(CTEST_MEMORYCHECK_COMMAND_OPTIONS "-v")
   smash_add_test(memcheck_box smash ${CTEST_MEMORYCHECK_COMMAND}
      ${CTEST_MEMORYCHECK_COMMAND_OPTIONS}
      ${PROJECT_BINARY_DIR}/smash -s 10 -m Box)
   smash_add_test(memcheck_collider smash ${CTEST_MEMORYCHECK_COMMAND}
      ${CTEST_MEMORYCHECK_COMMAND_OPTIONS}
      ${PROJECT_BINARY_DIR}/smash -s 50 -m Collider)
#  smash_add_test(memcheck_sphere smash ${CTEST_MEMORYCHECK_COMMAND}
#     ${CTEST_MEMORYCHECK_COMMAND_OPTIONS}
#     ${PROJECT_BINARY_DIR}/smash -s 50 -m Sphere)
else()
   message(STATUS "valgrind not found. Memcheck tests are disabled!")
endif()

if(HAVE_ADDRESS_SANITIZER)
   smash_add_test(asan_box      smash_asan smash_asan -s 500 -m Box)
   smash_add_test(asan_collider smash_asan smash_asan -s 500 -m Collider)
   #smash_add_test(asan_sphere   smash_asan smash_asan -s 500 -m Sphere)
endif()
