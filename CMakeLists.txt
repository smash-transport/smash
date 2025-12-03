########################################################
#
#    Copyright (c) 2012-2025
#      SMASH Team
#
#    BSD 3-clause license
#
########################################################

# Minimum cmake version this is tested on
cmake_minimum_required(VERSION 3.16)

# The name, version and language of our project
project(SMASH VERSION 3.3 LANGUAGES CXX)

# Fail if cmake is called in the source directory
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
    message(FATAL_ERROR "You don't want to configure in the source directory!")
endif()

# Restrict supported operating systems
if(NOT UNIX AND NOT APPLE OR CMAKE_SIZEOF_VOID_P LESS 8)
    option(FORCE_USE_ANY_SYSTEM "Force cmake to setup and compile on any OS/architecture." OFF)
    if(NOT FORCE_USE_ANY_SYSTEM)
        if(CMAKE_SIZEOF_VOID_P LESS 8)
            set(AUX_MSG "32-bit architecture")
        endif()
        if(NOT UNIX AND NOT APPLE)
            list(APPEND AUX_MSG "operating system")
            list(JOIN AUX_MSG "/" AUX_MSG)
        endif()
        message(FATAL_ERROR " \n" # See e.g. https://stackoverflow.com/a/51035045 for formatting
                            " Using SMASH on your ${AUX_MSG} is not officially supported.\n"
                            " Please, contact the SMASH development team (see README) if you would\n"
                            " like to collaborate with us to get the code run on your system.\n"
                            " At your own risk you can pass -DFORCE_USE_ANY_SYSTEM =TRUE to\n"
                            " cmake in order to setup and compile SMASH on your system.\n"
                            " Alternatively you can use the provided Docker containers (see README).\n"
        )
    endif()
endif()

# Restrict supported compilers
if(NOT CMAKE_CXX_COMPILER_ID MATCHES "^((Apple)?Clang|GNU)$")
    option(FORCE_USE_ANY_COMPILER "Force cmake to setup and compile with any compiler." OFF)
    if(NOT FORCE_USE_ANY_COMPILER)
        message(FATAL_ERROR " \n"
                            " Compiling SMASH with the specified compiler is not officially supported.\n"
                            " Please, contact the SMASH development team (see README) if you would\n"
                            " like to collaborate with us to get the codebase support your compiler.\n"
                            " At your own risk you can pass -DFORCE_USE_ANY_COMPILER=TRUE to cmake in\n"
                            " order to setup and try to compile SMASH with the chosen compiler.\n"
                            " Alternatively you can use the provided Docker containers (see README).\n"
        )
    endif()
endif()

# Before starting, tell cmake where to find our modules and include utilities
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(Utilities)

# The variable SMASH_VERSION is already set via the project() command, since we specify there a
# VERSION option. Here we check if .git information is available and increase verbosity
include(GetGitRevisionDescription)
git_describe(SMASH_VERSION_VERBOSE)
if(SMASH_VERSION_VERBOSE)
    set(SMASH_VERSION "${SMASH_VERSION_VERBOSE}")
else()
    set(SMASH_VERSION "${CMAKE_PROJECT_NAME}-${SMASH_VERSION}")
endif()

# Define installation top-level directory name inside which SMASH files are installed in sub-folders
string(TOLOWER "${SMASH_VERSION}" SMASH_INSTALLATION_SUBFOLDER)

# SMASH will be shipped as shared library and hence we want to build position independent code
include(CheckPIESupported)
check_pie_supported(LANGUAGES CXX)
if(NOT CMAKE_CXX_LINK_PIE_SUPPORTED)
    message(ATTENTION "PIE is not supported at link time for C++. "
                      "PIE link options will not be passed to linker.")
endif()
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Fail early if a too old officially supported compiler is used
if(CMAKE_CXX_COMPILER_ID MATCHES "^((Apple)?Clang|GNU)$")
    unset(MINIMUM_CXX_COMPILER_VERSION)
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        set(MINIMUM_CXX_COMPILER_VERSION "8")
    elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
        set(MINIMUM_CXX_COMPILER_VERSION "7")
    elseif(CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
        set(MINIMUM_CXX_COMPILER_VERSION "11")
    endif()
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS ${MINIMUM_CXX_COMPILER_VERSION})
        message(FATAL_ERROR " At least version ${MINIMUM_CXX_COMPILER_VERSION}"
                            " of ${CMAKE_CXX_COMPILER_ID} compiler is required.\n"
                            " Version in use is ${CMAKE_CXX_COMPILER_VERSION}")
    endif()
endif()

# Request a given C++ standard
set(CXX_SMASH_STANDARD "17")
if(NOT "${CMAKE_CXX_STANDARD}")
    set(CMAKE_CXX_STANDARD ${CXX_SMASH_STANDARD})
endif()
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# needed for clang-tidy
set(CMAKE_EXPORT_COMPILE_COMMANDS "ON")

# Set a default value for CMAKE_BUILD_TYPE, otherwise we get something which is none of the options.
if(NOT DEFINED CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "")
    set(CMAKE_BUILD_TYPE
        Release
        CACHE STRING
              "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel Profiling."
              FORCE)
elseif(NOT CMAKE_BUILD_TYPE MATCHES "^(Debug|Release|RelWithDebInfo|MinSizeRel|Profiling)$")
    message(FATAL_ERROR " \n"
                        " Invalid build configuration specified, CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}.\n"
                        " Valid values are: Debug, Release, RelWithDebInfo, MinSizeRel or Profiling.\n"
    )
elseif(CMAKE_BUILD_TYPE STREQUAL "Debug")
    message(ATTENTION "Consider adding '-O0' compiler flag for best debug information.")
endif()

# Retrieve architecture information about endianness
include(TestBigEndian)
test_big_endian(IS_BIG_ENDIAN)
if(IS_BIG_ENDIAN)
    message(STATUS "Big endian architecture detected.")
    option(FORCE_USE_ON_BIG_ENDIAN "Force cmake to setup and compile on big endian machine." OFF)
    if(NOT FORCE_USE_ON_BIG_ENDIAN)
        message(FATAL_ERROR " \n"
                            " Using SMASH on big endian machines is not officially supported.\n"
                            " Please, contact the SMASH development team (see README) if you would like\n"
                            " to collaborate with us to get the code run on big endian architectures.\n"
                            " At your own risk you can pass -DFORCE_USE_ON_BIG_ENDIAN=TRUE to cmake in\n"
                            " order to setup and compile SMASH on a big endian machine.\n"
                            " Alternatively you can use the provided Docker containers (see README).\n"
        )
    endif()
    add_definitions("-DBIG_ENDIAN_ARCHITECTURE")
else()
    message(STATUS "Little endian architecture detected.")
    add_definitions("-DLITTLE_ENDIAN_ARCHITECTURE")
endif()

# Let Clang users on Linux be able to use libc++
if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    option(CLANG_USE_LIBC++
           "If turned on clang will explicitly be asked to use libc++ (otherwise it uses the system default)"
           OFF)
    if(CLANG_USE_LIBC++)
        message(STATUS "Prepare compilation in order to use LLVM libc++ implementation")
        set(MESSAGE_QUIET ON)
        add_compiler_flag(-stdlib=libc++ CXX_FLAGS CMAKE_CXX_FLAGS CXX_RESULT _use_libcxx)
        unset(MESSAGE_QUIET)
        if(_use_libcxx AND "${CMAKE_SYSTEM_NAME}" STREQUAL "Linux")
            link_libraries(c++abi supc++)
        endif()
    endif()
endif()

# add 3rd-party libraries (before setting compiler flags etc)
include_directories(SYSTEM "${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/Cuba-4.2.2")
include_directories(SYSTEM "${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/einhard")
include_directories(SYSTEM "${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/yaml-cpp-0.8.0/include")
add_subdirectory(3rdparty)

# Set up compile options passed to all targets. From CMake documentation: "The flags in
# CMAKE_<LANG>_FLAGS will be passed to the compiler before those in the per-configuration
# CMAKE_<LANG>_FLAGS_<CONFIG> variant, and before flags added by the add_compile_options() or
# target_compile_options() commands."
add_compiler_flags_if_supported("-fno-math-errno" # Tell the compiler to ignore errno setting of
                                                  # math functions. This can help the compiler
                                                  # considerably in optimizing mathematical
                                                  # expressions.
                                "-march=native" # The '-march=native' flag is not supported e.g. on
                                                # Apple M1 machines
                                CXX_FLAGS CMAKE_CXX_FLAGS)

# CMake provides 4 build configurations (Debug, Release, RelWithDebInfo, MinSizeRel) together with
# sensible values for the corresponding CMAKE_(C|CXX)_FLAGS_<CONFIG> variables (cached and marked as
# advanced). Here we need to set those for Profiling build, only. Since the relevant variables are
# created and cached by CMake, we need to force their value in the cache if their value is false,
# e.g. empty.
#
# See https://cristianadam.eu/20190223/modifying-the-default-cmake-build-types/ for more infos.
if(NOT CMAKE_CXX_FLAGS_PROFILING)
    add_compiler_flags_if_supported("-O3" "-DNDEBUG" "-pg" CXX_FLAGS SUPPORTED_CXX_FLAGS_PROFILING)
endif()
if(NOT CMAKE_C_FLAGS_PROFILING)
    add_compiler_flags_if_supported("-O3" "-DNDEBUG" "-pg" C_FLAGS SUPPORTED_C_FLAGS_PROFILING)
endif()
if(NOT CMAKE_CXX_FLAGS_PROFILING)
    set(CMAKE_CXX_FLAGS_PROFILING "${SUPPORTED_CXX_FLAGS_PROFILING}"
        CACHE STRING "Flags used by the C++ compiler during profile builds." FORCE)
endif()
if(NOT CMAKE_C_FLAGS_PROFILING)
    set(CMAKE_C_FLAGS_PROFILING "${SUPPORTED_C_FLAGS_PROFILING}"
        CACHE STRING "Flags used by the C compiler during profile builds." FORCE)
endif()
if(NOT CMAKE_EXE_LINKER_FLAGS_PROFILING)
    set(CMAKE_EXE_LINKER_FLAGS_PROFILING "${CMAKE_EXE_LINKER_FLAGS_RELEASE} -pg"
        CACHE STRING "Flags used by the linker during profile builds." FORCE)
endif()
mark_as_advanced(CMAKE_CXX_FLAGS_PROFILING CMAKE_C_FLAGS_PROFILING CMAKE_EXE_LINKER_FLAGS_PROFILING)

# have binary in the build directory
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})

# enable standard CTest
include(CTest)

# subdirectories where the code is
add_subdirectory(src)
add_subdirectory(doc)

# Don’t make the install target depend on the all target (i.e. also all tests, if on)
set(CMAKE_SKIP_INSTALL_ALL_DEPENDENCY TRUE)

# Copy the default input files to the installation directory
install(DIRECTORY input/ DESTINATION "share/${SMASH_INSTALLATION_SUBFOLDER}/input_files")
