########################################################
#
#    Copyright (c) 2026
#      SMASH Team
#
#    BSD 3-clause license
#
########################################################

# Function to add a unit test to SMASH
function(smash_add_unit_test name)
    _smash_add_exe(${name})
    _smash_add_test(${name} "unit;code" ${name} ${name})
endfunction()

# Function to add an integration test to SMASH
function(smash_add_integration_test name)
    _smash_add_exe(${name})
    _smash_add_test(${name} "integration;code" ${name} ${name})
endfunction()

# Function to add a functional test to SMASH
function(smash_add_functional_test name)
    add_test(NAME functional_${name}
             COMMAND ${Python3_EXECUTABLE} "${PROJECT_SOURCE_DIR}/src/tests/functional/${name}.py"
                     --source "${PROJECT_SOURCE_DIR}" --binary "${PROJECT_BINARY_DIR}")
    set_tests_properties(functional_${name} PROPERTIES FIXTURES_REQUIRED fixture_compile_smash
                                                       LABELS "functional;physics")
    add_custom_target(run_functional_${name}
                      COMMAND ${Python3_EXECUTABLE}
                              "${PROJECT_SOURCE_DIR}/src/tests/functional/${name}.py" --source
                              "${PROJECT_SOURCE_DIR}" --binary "${PROJECT_BINARY_DIR}"
                      DEPENDS smash
                      COMMENT "Executing functional test ${name}"
                      VERBATIM)
endfunction()

# The following function add a physics "test" to run e.g. smash with certain arguments, choosing an
# output directory according to the test name
function(smash_add_run_test name depends)
    _smash_add_test(${name}
                    "run;physics"
                    ${depends}
                    ${ARGN}
                    -o
                    "${PROJECT_BINARY_DIR}/test_output/${name}"
                    -f)
endfunction()

# The following function sets up a Python virtual environment and installs the requirements for
# functional tests.  Since these are written in Python, here the interpreter is found, a virtual
# environment is setup, the interpreter of the virtual environment is looked for and the needed
# Python requirements are installed. If any step fails, functional tests are disabled. Note that it
# is not necessary to activate or deactivate the virtual environment because we run functional tests
# explicitly using the Python_EXECUTABLE which is now set to that of the virtual environment.
# ~~~
# Function "side effects" are the following:
#  - sets a variable passed as argument to TRUE if the setup was successful and to FALSE otherwise;
#  - sets the Python3_EXECUTABLE in the parent scope if the setup was successful, which is needed
#    to run the functional tests.
# ~~~
function(do_setup_for_functional_tests setup_was_successful)
    _smash_setup_python_venv(VENV_OK PYTHON_EXEC PYTHON_VENV_PATH)
    if(VENV_OK)
        _smash_install_python_requirements(${PYTHON_EXEC} ${PYTHON_VENV_PATH} REQUIREMENTS_OK)
    else()
        set(setup_was_successful FALSE PARENT_SCOPE)
    endif()
    if(REQUIREMENTS_OK)
        # Since all functional tests need smash executable to be compiled before running, we want
        # ctest to first compile SMASH and then run the functional tests. In order to do so we use
        # tests fixtures setting up a "test" to compile SMASH and using this as setup fixture. See
        # CMake doc https://cmake.org/cmake/help/latest/prop_test/FIXTURES_REQUIRED.html for further
        # information. Please note that we do not want to always compile SMASH before each test and
        # hence we wrap the SMASH compilation into a small bash test to skip it if the smash
        # executable exists.
        _smash_setup_functional_tests_fixture()
        set(Python3_EXECUTABLE ${PYTHON_EXEC} PARENT_SCOPE)
        set(setup_was_successful TRUE PARENT_SCOPE)
    else()
        set(setup_was_successful FALSE PARENT_SCOPE)
    endif()
endfunction()

# =================================================================================================
# The following part contains the implementation details of the functionality called from the
# functions above and should not be used directly!
# =================================================================================================

# The following function adds an executable for a test, linking it to the test library and taking
# care of some compiler options to be adjusted.
function(_smash_add_exe name)
    add_executable(${name} EXCLUDE_FROM_ALL ${name}.cc)
    add_dependencies(tests ${name})
    target_link_libraries(${name} PRIVATE smash_testlib)
    set_target_properties(${name} PROPERTIES CXX_EXTENSIONS OFF)
    target_compile_definitions(${name}
                               PRIVATE SMASH_TEST_OUTPUT_PATH="${PROJECT_BINARY_DIR}/test_output/${name}"
    )
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        # Disable variable tracking in debug mode (never used so far and gcc gives annoying notes
        # when failing to track them)
        target_compile_options(${name} PRIVATE "-fno-var-tracking-assignments")
    endif()
    if(CMAKE_CXX_COMPILER_ID MATCHES "^((Apple)?Clang|GNU)$")
        # It turned out that the -ffp-contract option of GNU and LLVM compilers (whose default is
        # different in the two compilers and among versions) is affecting reproducibility of results
        # when e.g. optimizing for different architectures (e.g. -march=native VS -march=x86-64).
        # Therefore, we decided to switch this feature off for tests to ensure that tests behave in
        # the same way on all machines.
        target_compile_options(${name} PRIVATE "-ffp-contract=off")
    endif()
endfunction()

# The following function adds a test associated to the <name> executable and creates a run_<name>
# target with the correct dependencies. This enables running the test e.g. using make. It is meant
# to be an implementation detail and the code adding test should not use this function.
function(_smash_add_test name label depends)
    add_test(NAME ${name} COMMAND ${ARGN})
    set_tests_properties(${name} PROPERTIES LABELS "${label}")
    if(APPLE)
        # Known false positive container overflow error of AddressSanitizer on Mac.
        # https://github.com/google/sanitizers/wiki/AddressSanitizerContainerOverflow
        set_tests_properties(${name} PROPERTIES ENVIRONMENT
                                                "ASAN_OPTIONS=detect_container_overflow=0")
    endif()
    add_custom_target(run_${name}
                      COMMAND ${ARGN}
                      DEPENDS ${depends}
                      COMMENT "Executing test ${name}"
                      VERBATIM)
endfunction()

# The following function sets up a Python virtual environment which is needed for functional tests
function(_smash_setup_python_venv out_success out_python_exec out_venv_path)
    message(STATUS "Looking for Python3")
    # Note that the functional tests require Pandas 2.0 which in turn requires at least Python 3.8.
    # However, Python 3.12 (and above) is excluded because it contains setuptools>=82 and with it
    # pkg_resources has been removed which causes problems installing Pandas v2.x.
    find_package(Python3 3.8...<3.12 QUIET COMPONENTS Interpreter Development)
    if(NOT Python3_FOUND)
        message(ATTENTION "Python3 not found. Functional tests disabled.")
        set(${out_success} OFF PARENT_SCOPE)
        return()
    endif()
    message(STATUS "Found Python3: ${Python3_EXECUTABLE}")
    set(PYTHON_VENV_PATH "${PROJECT_BINARY_DIR}/python_venv")
    execute_process(COMMAND ${Python3_EXECUTABLE} -m venv "${PYTHON_VENV_PATH}"
                    RESULT_VARIABLE venv_result)
    if(venv_result AND NOT venv_result EQUAL 0)
        message(ATTENTION "Setting up Python3 venv failed. Functional tests disabled.")
        set(${out_success} OFF PARENT_SCOPE)
        return()
    endif()
    # Look again for Python executable, this time to locate that of the venv, which is needed to
    # install the requirements and run the functional tests (taken from CMake discourse)
    set(SYSTEM_PYTHON_EXECUTABLE ${Python3_EXECUTABLE})
    set(ENV{VIRTUAL_ENV} "${PYTHON_VENV_PATH}") # Update the environment with VIRTUAL_ENV variable
                                                # (mimic the activate script)
    set(Python3_FIND_VIRTUALENV FIRST) # Change the context of the search
    unset(Python3_EXECUTABLE) # Needed as it is also an input variable for find_package
    find_package(Python3 QUIET COMPONENTS Interpreter Development)
    if(SYSTEM_PYTHON_EXECUTABLE STREQUAL Python3_EXECUTABLE)
        message(FATAL_ERROR " \n" " Python3 venv executable is the same as the system one.\n"
                            " This should not happen and needs to be investigated.\n")
    else()
        message(STATUS "Python3 exec for functional tests: ${Python3_EXECUTABLE}")
        message(STATUS "Python3 venv for functional tests: ${PYTHON_VENV_PATH}")
        set(${out_success} ON PARENT_SCOPE)
        set(${out_python_exec} ${Python3_EXECUTABLE} PARENT_SCOPE)
        set(${out_venv_path} ${PYTHON_VENV_PATH} PARENT_SCOPE)
    endif()
endfunction()

# The following function installs the Python requirements for functional tests into the venv
function(_smash_install_python_requirements python_exec venv_path out_success)
    message(STATUS "Installing Python3 requirements for functional tests into venv")
    execute_process(COMMAND ${python_exec} -m pip install --upgrade pip -r
                            "${PROJECT_SOURCE_DIR}/src/tests/functional/requirements.txt"
                    RESULT_VARIABLE pip_install_result
                    OUTPUT_FILE "${venv_path}/pip_install_requirements.log")
    if(pip_install_result AND NOT pip_install_result EQUAL 0)
        message(ATTENTION "Installing Python3 requirements failed. Functional tests disabled.")
        set(${out_success} OFF PARENT_SCOPE)
    else()
        message(STATUS "Installing Python3 requirements for functional tests into venv - done")
        set(${out_success} ON PARENT_SCOPE)
    endif()
endfunction()

# The following function adds a test to compile the SMASH executable, which is needed for functional
# tests, and sets up a test fixture to ensure that it is run before any functional test
function(_smash_setup_functional_tests_fixture)
    if(NOT TEST compile_smash_executable)
        set(BASH_CODE_FOR_TEST
            "[[ ! -f '${CMAKE_BINARY_DIR}/smash' ]] && \
            ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --config ${CMAKE_BUILD_TYPE} \
                             --target smash  || exit 0")
        add_test(NAME compile_smash_executable COMMAND bash -c ${BASH_CODE_FOR_TEST})
        set_tests_properties(compile_smash_executable PROPERTIES FIXTURES_SETUP
                                                                 fixture_compile_smash)
    endif()
endfunction()
