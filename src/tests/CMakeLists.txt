# test sources
add_executable(fourvector EXCLUDE_FROM_ALL fourvector.cc)
add_executable(lorentzboost EXCLUDE_FROM_ALL lorentzboost.cc ../fourvector.cc)
add_executable(angles EXCLUDE_FROM_ALL angles.cc)
add_executable(angles_distribution EXCLUDE_FROM_ALL angles_distribution.cc)
add_executable(particles EXCLUDE_FROM_ALL particles.cc ../particles.cc ../fourvector.cc ../distributions.cc ../outputroutines.cc ../time.cc)
add_executable(clebschgordan EXCLUDE_FROM_ALL clebschgordan.cc ../clebschgordan.cc)
add_dependencies(check fourvector lorentzboost particles clebschgordan angles
angles_distribution smash)

file(COPY ../config_general.txt DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ../particles.txt DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ../decaymodes.txt DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

target_link_libraries(particles ${SMASH_LIBRARIES})


# tests:

# test source code to be clean, requires cpplint in the path
file(GLOB source "../*.cc" "../include/*.h" "*.cc")
add_test(cclint cpplint.py ${source})

# some fourvector operations tests
add_test(fourvector ${PROJECT_BINARY_DIR}/fourvector)

# Lorentz boost tests
add_test(lorentzboost ${PROJECT_BINARY_DIR}/lorentzboost)

# some angle operations tests
add_test(angles ${PROJECT_BINARY_DIR}/angles)

# some particles operations tests
add_test(particles ${PROJECT_BINARY_DIR}/particles)

# clebsch-gordan coefficient tests
add_test(clebschgordan ${PROJECT_BINARY_DIR}/clebschgordan)

# verify that the binary has a cli help
add_test(mash_help ${PROJECT_BINARY_DIR}/smash -h)

# verify that the binary runs
add_test(mash_run ${PROJECT_BINARY_DIR}/smash)

# verify that the binary can run a box calculation
set(CTEST_SMASH_OPTIONS -m Box)
add_test(box_run ${PROJECT_BINARY_DIR}/smash ${CTEST_SMASH_OPTIONS})
# verify that the binary can do Collider
set(CTEST_SMASH_OPTIONS -m Collider)
add_test(collider_run ${PROJECT_BINARY_DIR}/smash ${CTEST_SMASH_OPTIONS})
# verify that the binary can do sphere
#set(CTEST_SMASH_OPTIONS -m Sphere)
#add_test(sphere_run ${PROJECT_BINARY_DIR}/smash ${CTEST_MASH_OPTIONS})

# verify that binary runs certain number of steps
set(CTEST_SMASH_OPTIONS -s 500 -m Box)
add_test(box_steps ${PROJECT_BINARY_DIR}/smash ${CTEST_SMASH_OPTIONS})
set(CTEST_SMASH_OPTIONS -s 500 -m Collider)
add_test(collider_steps ${PROJECT_BINARY_DIR}/smash ${CTEST_SMASH_OPTIONS})
#set(CTEST_MASH_COMMAND_OPTIONS -S 500 -m Sphere)
#add_test(sphere_steps ${PROJECT_BINARY_DIR}/smash ${CTEST_MASH_COMMAND_OPTIONS})

# verify that default run has no mem leaks
set(CTEST_MEMORYCHECK_COMMAND valgrind )
set(CTEST_MEMORYCHECK_COMMAND_OPTIONS "-v")
set(CTEST_SMASH_OPTIONS -s 500 -m Box)
add_test(memcheck_box ${CTEST_MEMORYCHECK_COMMAND} ${CTEST_MEMORYCHECK_COMMAND_OPTIONS} ${PROJECT_BINARY_DIR}/smash ${CTEST_SMASH_OPTIONS})
set(CTEST_SMASH_OPTIONS -s 500 -m Collider)
add_test(memcheck_collider ${CTEST_MEMORYCHECK_COMMAND} ${CTEST_MEMORYCHECK_COMMAND_OPTIONS} ${PROJECT_BINARY_DIR}/smash ${CTEST_SMASH_OPTIONS})
#set(CTEST_SMASH_OPTIONS -S 500 -m Sphere)
#add_test(memcheck_sphere ${CTEST_MEMORYCHECK_COMMAND} ${CTEST_MEMORYCHECK_COMMAND_OPTIONS} ${PROJECT_BINARY_DIR}/smash ${CTEST_SMASH_OPTIONS})
