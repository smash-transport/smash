# test sources
add_executable(fourvector fourvector.cc)
add_executable(lorentzboost lorentzboost.cc ../fourvector.cc)
add_executable(angles angles.cc)
add_executable(angles_distribution angles_distribution.cc)
add_executable(particles particles.cc ../particles.cc ../fourvector.cc ../distributions.cc ../outputroutines.cc ../time.cc)
add_executable(clebschgordan clebschgordan.cc ../clebschgordan.cc)

file(COPY ../config_general.txt DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ../particles.txt DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ../decaymodes.txt DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

target_link_libraries(particles ${SMASH_LIBRARIES})


# tests:
macro(smash_add_test name depends)
   add_test(NAME ${name} COMMAND ${ARGN})
   add_custom_target(run_${name}
      COMMAND ${ARGN}
      DEPENDS ${depends}
      COMMENT "Executing test ${name}"
      VERBATIM
      )
endmacro()

# test source code to be clean, requires cpplint in the path
file(GLOB source "../*.cc" "../include/*.h" "*.cc")
add_test(cclint cpplint.py ${source})

# some fourvector operations tests
smash_add_test(fourvector fourvector fourvector)

# Lorentz boost tests
smash_add_test(lorentzboost lorentzboost lorentzboost)

# some angle operations tests
smash_add_test(angles angles angles)

# some particles operations tests
smash_add_test(particles particles particles)

# clebsch-gordan coefficient tests
smash_add_test(clebschgordan clebschgordan clebschgordan)

# verify that the binary has a cli help
smash_add_test(mash_help smash smash -h)

# verify that the binary runs
smash_add_test(mash_run smash smash)

# verify that the binary can run a box calculation
set(CTEST_SMASH_OPTIONS -m Box)
smash_add_test(box_run smash smash ${CTEST_SMASH_OPTIONS})
# verify that the binary can do Collider
set(CTEST_SMASH_OPTIONS -m Collider)
smash_add_test(collider_run smash smash ${CTEST_SMASH_OPTIONS})
# verify that the binary can do sphere
#set(CTEST_SMASH_OPTIONS -m Sphere)
#smash_add_test(sphere_run smash smash ${CTEST_MASH_OPTIONS})

# verify that binary runs certain number of steps
set(CTEST_SMASH_OPTIONS -s 500 -m Box)
smash_add_test(box_steps smash smash ${CTEST_SMASH_OPTIONS})
set(CTEST_SMASH_OPTIONS -s 500 -m Collider)
smash_add_test(collider_steps smash smash ${CTEST_SMASH_OPTIONS})
#set(CTEST_MASH_COMMAND_OPTIONS -S 500 -m Sphere)
#smash_add_test(sphere_steps smash smash ${CTEST_MASH_COMMAND_OPTIONS})

# verify that default run has no mem leaks
set(CTEST_MEMORYCHECK_COMMAND valgrind )
set(CTEST_MEMORYCHECK_COMMAND_OPTIONS "-v")
set(CTEST_SMASH_OPTIONS -s 500 -m Box)
smash_add_test(memcheck_box smash ${CTEST_MEMORYCHECK_COMMAND} ${CTEST_MEMORYCHECK_COMMAND_OPTIONS} ${PROJECT_BINARY_DIR}/smash ${CTEST_SMASH_OPTIONS})
set(CTEST_SMASH_OPTIONS -s 500 -m Collider)
smash_add_test(memcheck_collider smash ${CTEST_MEMORYCHECK_COMMAND} ${CTEST_MEMORYCHECK_COMMAND_OPTIONS} ${PROJECT_BINARY_DIR}/smash ${CTEST_SMASH_OPTIONS})
#set(CTEST_SMASH_OPTIONS -S 500 -m Sphere)
#smash_add_test(memcheck_sphere smash ${CTEST_MEMORYCHECK_COMMAND} ${CTEST_MEMORYCHECK_COMMAND_OPTIONS} ${PROJECT_BINARY_DIR}/smash ${CTEST_SMASH_OPTIONS})

if(HAVE_ADDRESS_SANITIZER)
   smash_add_test(asan_box      smash_asan smash_asan -s 500 -m Box)
   smash_add_test(asan_collider smash_asan smash_asan -s 500 -m Collider)
   #smash_add_test(asan_sphere   smash_asan smash_asan -s 500 -m Sphere)
endif()
