add_executable(fourvector EXCLUDE_FROM_ALL fourvector.cc)
add_executable(particles EXCLUDE_FROM_ALL particles.cc ../particles.cc ../fourvector.cc)
add_dependencies(check fourvector particles)

# tests:

# test source code to be clean, requires cpplint in the path
file(GLOB source "../*.cc" "../include/*.h" "*.cc")
add_test(cclint cpplint.py ${source})

# some fourvector operations tests
add_test(fourvector ${PROJECT_BINARY_DIR}/fourvector)

# some particles operations tests
add_test(particles ${PROJECT_BINARY_DIR}/particles)

# verify that the binary has a cli help
add_test(mash_help ${PROJECT_BINARY_DIR}/mash -h)

# verify that the binary runs
add_test(mash_run ${PROJECT_BINARY_DIR}/mash)

# verify that binary runs certain number of steps
set(CTEST_MASH_COMMAND_OPTIONS "-s 500")
add_test(mash_steps ${PROJECT_BINARY_DIR}/mash ${CTEST_MASH_COMMAND_OPTIONS})

# verify that default run has no mem leaks
set(CTEST_MEMORYCHECK_COMMAND valgrind )
set(CTEST_MEMORYCHECK_COMMAND_OPTIONS "-v")
add_test(memcheck ${CTEST_MEMORYCHECK_COMMAND} ${CTEST_MEMORYCHECK_COMMAND_OPTIONS} ${PROJECT_BINARY_DIR}/mash ${CTEST_MASH_COMMAND_OPTIONS})
