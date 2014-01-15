# test sources
add_executable(fourvector EXCLUDE_FROM_ALL fourvector.cc)
add_executable(particles EXCLUDE_FROM_ALL particles.cc ../particles.cc ../fourvector.cc ../distributions.cc ../outputroutines.cc ../time.cc)
add_executable(clebschgordan EXCLUDE_FROM_ALL clebschgordan.cc ../clebschgordan.cc)
add_dependencies(check fourvector particles clebschgordan)

# link flags: rt may be needed for old glibc
if(HAVE_CLOCK_GETTIME)
  target_link_libraries(particles gsl gslcblas rt)
else(HAVE_CLOCK_GETTIME)
  target_link_libraries(particles gsl gslcblas)
endif(HAVE_CLOCK_GETTIME)


# tests:

# test source code to be clean, requires cpplint in the path
file(GLOB source "../*.cc" "../include/*.h" "*.cc")
add_test(cclint cpplint.py ${source})

# some fourvector operations tests
add_test(fourvector ${PROJECT_BINARY_DIR}/fourvector)

# some particles operations tests
add_test(particles ${PROJECT_BINARY_DIR}/particles)

# clebsch-gordan coefficient tests
add_test(clebschgordan ${PROJECT_BINARY_DIR}/clebschgordan)

# verify that the binary has a cli help
add_test(mash_help ${PROJECT_BINARY_DIR}/smash -h)

# verify that the binary runs
add_test(mash_run ${PROJECT_BINARY_DIR}/smash)

# verify that the binary can do sphere
# verify that binary runs certain number of steps
set(CTEST_MASH_COMMAND_OPTIONS "-m 2")
add_test(sphere_run ${PROJECT_BINARY_DIR}/smash ${CTEST_MASH_COMMAND_OPTIONS})

# verify that binary runs certain number of steps
set(CTEST_MASH_COMMAND_OPTIONS "-S 500")
add_test(mash_steps ${PROJECT_BINARY_DIR}/smash ${CTEST_MASH_COMMAND_OPTIONS})
set(CTEST_MASH_COMMAND_OPTIONS "-S 500 -m 2")
add_test(sphere_steps ${PROJECT_BINARY_DIR}/smash ${CTEST_MASH_COMMAND_OPTIONS})

# verify that default run has no mem leaks
set(CTEST_MEMORYCHECK_COMMAND valgrind )
set(CTEST_MEMORYCHECK_COMMAND_OPTIONS "-v")
set(CTEST_MASH_COMMAND_OPTIONS "-S 500")
add_test(memcheck ${CTEST_MEMORYCHECK_COMMAND} ${CTEST_MEMORYCHECK_COMMAND_OPTIONS} ${PROJECT_BINARY_DIR}/smash ${CTEST_MASH_COMMAND_OPTIONS})
set(CTEST_MASH_COMMAND_OPTIONS "-S 500 -m 2")
add_test(memcheck_sphere ${CTEST_MEMORYCHECK_COMMAND} ${CTEST_MEMORYCHECK_COMMAND_OPTIONS} ${PROJECT_BINARY_DIR}/smash ${CTEST_MASH_COMMAND_OPTIONS})
