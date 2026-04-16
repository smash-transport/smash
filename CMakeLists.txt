########################################################
#
#    Copyright (c) 2012-2026
#      SMASH Team
#
#    BSD 3-clause license
#
########################################################

# Minimum cmake version this is tested on
cmake_minimum_required(VERSION 3.28)

# The name, version and language of our project
project(SMASH VERSION 3.3 LANGUAGES CXX)
set(SMASH_VERSION "${SMASH_VERSION}-next")

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
                            " At your own risk you can pass -DFORCE_USE_ANY_SYSTEM=TRUE to\n"
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
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS MINIMUM_CXX_COMPILER_VERSION)
        message(FATAL_ERROR " At least version ${MINIMUM_CXX_COMPILER_VERSION}"
                            " of ${CMAKE_CXX_COMPILER_ID} compiler is required.\n"
                            " Version in use is ${CMAKE_CXX_COMPILER_VERSION}")
    endif()
endif()

# Before starting, tell cmake where to find our modules and include utilities
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
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
set(CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}/${SMASH_INSTALLATION_SUBFOLDER}"
    CACHE PATH "Install prefix for the present SMASH version" FORCE)

# SMASH will be shipped as shared library and hence we want to build position independent code
include(CheckPIESupported)
check_pie_supported(LANGUAGES CXX)
if(NOT CMAKE_CXX_LINK_PIE_SUPPORTED)
    message(ATTENTION "PIE is not supported at link time for C++. "
                      "PIE link options will not be passed to linker.")
endif()
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Generate compile_commands.json for better integration with tools like clang-tidy. Adding it to the
# cache makes it visible to IDEs or tools consuming the build.
set(CMAKE_EXPORT_COMPILE_COMMANDS ON CACHE BOOL "Generate compile_commands.json" FORCE)

# Set a default value for CMAKE_BUILD_TYPE, otherwise we get something which is none of the options.
if(NOT DEFINED CMAKE_BUILD_TYPE OR NOT CMAKE_BUILD_TYPE)
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

# Define two interface libraries to hold usage and build requirements of the project. Build
# requirements are those internal to SMASH e.g. compiler flags that we want to add to all targets,
# while usage requirements are all properties which are also needed by SMASH clients.
add_library(smash_project_build_requirements INTERFACE)
add_library(smash_project_usage_requirements INTERFACE)

# Request a given C++ standard. Note that CMake will error if the compiler does not support the
# requested standard, without any need to set the global CMAKE_CXX_STANDARD_REQUIRED variable. We
# would like to also disable here compiler-specific extensions, but properties like CXX_EXTENSIONS
# do not propagate through INTERFACE libraries and hence we need to set it on the targets directly.
target_compile_features(smash_project_usage_requirements INTERFACE cxx_std_17)

# Retrieve architecture information about endianness. Note that this compile definition is ATM a
# usage requirement as the PDG code header file relies on it.
if(CMAKE_CXX_BYTE_ORDER STREQUAL "BIG_ENDIAN")
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
    target_compile_definitions(smash_project_usage_requirements INTERFACE BIG_ENDIAN_ARCHITECTURE)
else()
    message(STATUS "Little endian architecture detected.")
    target_compile_definitions(smash_project_usage_requirements
                               INTERFACE LITTLE_ENDIAN_ARCHITECTURE)
endif()

# Let Clang users on Linux be able to use libc++. Note that on Apple machines libc++ is the default
# and hence this is not needed. Warn the user, if they try to enable this option.
option(CLANG_USE_LIBC++
       "If turned on clang will explicitly be asked to use libc++ (otherwise it uses the system default)"
       OFF)
if(CLANG_USE_LIBC++)
    if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
        message(STATUS "Using LLVM libc++ implementation for ${CMAKE_CXX_COMPILER_ID} compiler")
        target_add_compiler_flag_if_supported(TARGETS smash_project_build_requirements
                                              LANGUAGE CXX
                                              FLAGS -stdlib=libc++
                                              SCOPE INTERFACE RESULT _libcxx_supported)
        list(GET _libcxx_supported 0 _use_libcxx)
        if(_use_libcxx AND "${CMAKE_SYSTEM_NAME}" STREQUAL "Linux")
            message(STATUS "Linking libc++ support libraries for Linux Clang")
            target_link_libraries(smash_project_build_requirements INTERFACE c++abi supc++)
        endif()
    else()
        set(_message_text
            " CLANG_USE_LIBC++ option set to ON, but the compiler is not Clang; ignoring it.\n")
        if(CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
            string(APPEND
                   _message_text
                   " Note that AppleClang on macOS uses libc++ as the default C++ standard library,\n"
                   " so this option is not needed there.\n")
        endif()
        message(WARNING " \n" "${_message_text}")
    endif()
endif()

# Set up compile options that are considered build requirements for SMASH. These flags are attached
# to smash_project_build_requirements INTERFACE, so they propagate to all targets that link to them.
# Since CMake applies flags in this order (as documented):
# ~~~
# (1) CMAKE_<LANG>_FLAGS
# (2) CMAKE_<LANG>_FLAGS_<CONFIG>
# (3) add_compile_options()
# (4) target_compile_options()
# ~~~
# we want to provide a way for users to change the architecture tuning, which is by default set to
# native. Given the order above, any attempt to set architecture tuning via CMAKE_CXX_FLAGS will be
# overridden by the value set here and hence we want to warn users if they try to do so. Furthermore
# we want to provide a cache variable for users to properly set the architecture tuning, which is
# then passed to the compiler via target_compile_options().
set(SMASH_ARCH_FLAG
    "native"
    CACHE STRING
          "Value passed to -march compiler flag (e.g. native, x86-64). Empty string disables architecture tuning."
)
# Warn the user if they try to set architecture tuning via CMAKE_CXX_FLAGS, which is not the correct
# way and will be overridden by SMASH_ARCH_FLAG.
if(CMAKE_CXX_FLAGS MATCHES "-march")
    message(WARNING " \n"
                    " Detected '-march=' inside CMAKE_CXX_FLAGS. SMASH controls\n"
                    " architecture tuning via SMASH_ARCH_FLAG. This flag provided\n"
                    " via CMAKE_CXX_FLAGS will be overridden. Please use:\n"
                    "   -DSMASH_ARCH_FLAG=<value>\n"
                    " instead of modifying CMAKE_CXX_FLAGS.")
endif()
# Finally, add the architecture tuning flag if supported by the compiler and inform the user about
# the current value of architecture tuning and how to change it.
set(_architecture_tuning_message "Architecture tuning controlled by SMASH_ARCH_FLAG (current:")
if(SMASH_ARCH_FLAG)
    message(ATTENTION "${_architecture_tuning_message} -march=${SMASH_ARCH_FLAG})")
    target_add_compiler_flag_if_supported(TARGETS smash_project_build_requirements
                                          LANGUAGE CXX
                                          FLAGS -march=${SMASH_ARCH_FLAG}
                                          SCOPE INTERFACE)
else()
    message(STATUS "${_architecture_tuning_message} disabled)")
endif()
# Tell the compiler to ignore errno setting of math functions. This can help the compiler
# considerably in optimizing mathematical expressions.
target_add_compiler_flag_if_supported(TARGETS smash_project_build_requirements
                                      LANGUAGE CXX
                                      FLAGS -fno-math-errno
                                      SCOPE INTERFACE)

# CMake provides 4 build configurations (Debug, Release, RelWithDebInfo, MinSizeRel) together with
# sensible values for the corresponding CMAKE_(C|CXX)_FLAGS_<CONFIG> variables (cached and marked as
# advanced). Here we need to set those for Profiling build, only. We only set a default if the cache
# variable is not defined, allowing users to override it (even with an empty value). However, we
# warn users if they set an empty value, since that might lead to non-instrumented builds which are
# not useful for profiling.
#
# See https://cristianadam.eu/20190223/modifying-the-default-cmake-build-types/ for more infos.
#
# Note that here we do not test whether the compiler supports the provided flags, because these are
# universal for compilers that we support. Do not set up Profiling build flags for other compilers
# and warn user about them having to set up those flags manually.
if(CMAKE_CXX_COMPILER_ID MATCHES "^((Apple)?Clang|GNU)$")
    if(NOT DEFINED CMAKE_CXX_FLAGS_PROFILING)
        set(CMAKE_CXX_FLAGS_PROFILING "-O3 -DNDEBUG -pg"
            CACHE STRING "Flags used by the C++ compiler during Profiling builds.")
    elseif(CMAKE_CXX_FLAGS_PROFILING STREQUAL "")
        message(WARNING "CMAKE_CXX_FLAGS_PROFILING is empty. "
                        "Profiling builds may not generate useful profiling data.")
    endif()
    if(NOT DEFINED CMAKE_C_FLAGS_PROFILING)
        set(CMAKE_C_FLAGS_PROFILING "-O3 -DNDEBUG -pg"
            CACHE STRING "Flags used by the C compiler during Profiling builds.")
    elseif(CMAKE_C_FLAGS_PROFILING STREQUAL "")
        message(WARNING "CMAKE_C_FLAGS_PROFILING is empty. "
                        "Profiling builds may not generate useful profiling data.")
    endif()
    if(NOT DEFINED CMAKE_EXE_LINKER_FLAGS_PROFILING)
        set(CMAKE_EXE_LINKER_FLAGS_PROFILING "${CMAKE_EXE_LINKER_FLAGS_RELEASE} -pg"
            CACHE STRING "Flags used by the linker during Profiling builds.")
    elseif(CMAKE_EXE_LINKER_FLAGS_PROFILING STREQUAL "")
        message(WARNING "CMAKE_EXE_LINKER_FLAGS_PROFILING is empty. "
                        "Profiling builds may not generate useful profiling data.")
    endif()
    mark_as_advanced(CMAKE_CXX_FLAGS_PROFILING CMAKE_C_FLAGS_PROFILING
                     CMAKE_EXE_LINKER_FLAGS_PROFILING)
else()
    message(WARNING " \n"
                    " Profiling build flags are not set for compiler ${CMAKE_CXX_COMPILER_ID}.\n"
                    " You have to set them manually, if you need a Profiling build.")
endif()
mark_as_advanced(CMAKE_CXX_FLAGS_PROFILING CMAKE_C_FLAGS_PROFILING CMAKE_EXE_LINKER_FLAGS_PROFILING)

# By default, produce binary in the build directory. Note that the target property
# RUNTIME_OUTPUT_DIRECTORY supersedes this variable for a target if it is set.
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})

# Enable standard CTest
include(CTest)

# Add subdirectories where the code is. Note that SMASH requirements are propagated through targets
# and no global change should have occurred above. Therefore, adding the vendored third-party
# libraries here should not cause any issue, since they will not be affected by SMASH specific build
# or usage requirements. On the contrary, we can use there the SMASH_ARCH_FLAG variable to set the
# architecture tuning for the compilation of the vendored libraries as well, which is necessary.
add_subdirectory(3rdparty)
add_subdirectory(src)
add_subdirectory(doc)

# Start setting up installation of SMASH components as well as CMake configuration version file for
# find_package. Note that we install the non-standard modules to find dependencies as well, since
# they are needed by users of SMASH as a library. All CMake related files are installed in
# share/cmake/SMASH, which is a common practice for CMake configuration files. We also install the
# default input files, since they might be needed by users running SMASH. The SMASHConfig.cmake
# configuration file will be configured and installed in the src/CMakeLists.txt file, which is the
# place where all information will be available. For the moment we use AnyNewerVersion compatibility
# for the version file, since we do not have a stable API yet and we want to allow users to use
# newer versions of SMASH without having to update their CMake configuration files. Once we have a
# stable API, we can change this to SameMajorVersion or stricter options, depending on the level of
# compatibility we want to guarantee.
include(CMakePackageConfigHelpers)
write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/SMASHConfigVersion.cmake
                                 VERSION ${PROJECT_VERSION} COMPATIBILITY AnyNewerVersion)
install(EXPORT SMASH-targets-core
        FILE SMASH-targets-core.cmake
        NAMESPACE SMASH::
        DESTINATION share/cmake/SMASH)
install(FILES ${CMAKE_SOURCE_DIR}/cmake/FindPythia.cmake ${CMAKE_SOURCE_DIR}/cmake/FindEigen3.cmake
              ${CMAKE_CURRENT_BINARY_DIR}/SMASHConfigVersion.cmake DESTINATION share/cmake/SMASH)

# Copy the default input files to the installation directory
install(DIRECTORY input/ DESTINATION "share/input_files")
