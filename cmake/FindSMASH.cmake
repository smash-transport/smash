########################################################
#
#    Copyright (c) 2018-2020,2022-2025
#      SMASH Team
#
#    BSD 3-clause license
#
########################################################

# cmake-format: off
#===============================================================================
# - Try to find SMASH installation
#
# Once done this will define:
#
# SMASH_LIBRARIES         smash libraries
# SMASH_INCLUDE_DIR       directory of smash includes
# SMASH_INPUT_FILES_DIR   directory containing example input files
# SMASH_FOUND             true is smash found
#
# This module gives more informative information about the outcome if run
# using the --log-level=VERBOSE CMake command line option.
#
# The environment variable SMASH_INSTALL_DIR must be properly set to succeed
# and, in particular, refer to where SMASH has been intalled via 'make install'.
# For instance, use
#   export SMASH_INSTALL_DIR="${HOME}/.local/"
# if SMASH has been installed via
#   cmake -DCMAKE_INSTALL_PREFIX="${HOME}/.local" ..
#   make install
#
# Optionally, a SMASH_VERSION environment variable can be defined to specify
# which version should be found. This will be used as suffix after 'smash-' to
# locate files in the SMASH_INSTALL_DIR sub-folders.
#
# DEPRECATED APPROACH:
# The environment variable SMASH_DIR must be set properly to succeed, e.g.:
#   export SMASH_DIR=~/Work/SMASH/smash
#
# Optionally, a SMASH_BUILD_DIR environment variable can be defined to specify
# the path to the folder where SMASH has been built. If this is not set, it is
# assumed that SMASH has been built in "${SMASH_DIR}/build".
#===============================================================================
# cmake-format: on

# At the moment a macro about the system endianness is needed from within SMASH
include(TestBigEndian)
test_big_endian(IS_BIG_ENDIAN)
if(IS_BIG_ENDIAN)
    message(STATUS "Big endian architecture detected.")
    add_definitions("-DBIG_ENDIAN_ARCHITECTURE")
else()
    message(STATUS "Little endian architecture detected.")
    add_definitions("-DLITTLE_ENDIAN_ARCHITECTURE")
endif()

message(STATUS "Looking for SMASH")

if(NOT DEFINED ENV{SMASH_INSTALL_DIR})
    if(DEFINED ENV{SMASH_DIR})
        message(WARNING " \n"
                        " Usage of the SMASH_DIR environment variable containing the path to\n"
                        " the SMASH codebase is deprecated in the FindSMASH.cmake module.\n"
                        " Support for it might be removed in the future. Prefer installing\n"
                        " SMASH and then using SMASH_INSTALL_DIR in other projects.\n")
        if(NOT DEFINED ENV{SMASH_BUILD_DIR})
            set(SMASH_BUILD_DIR "$ENV{SMASH_DIR}/build")
        else()
            set(SMASH_BUILD_DIR "$ENV{SMASH_BUILD_DIR}")
        endif()
    else()
        message(FATAL_ERROR " \n"
                            " Environment variable SMASH_INSTALL_DIR containing the path to the\n"
                            " SMASH installation found undefined but needed by FindSMASH.cmake module.\n"
        )
    endif()
else()
    if(DEFINED ENV{SMASH_DIR} OR DEFINED ENV{SMASH_BUILD_DIR})
        message(WARNING " \n" " SMASH_DIR and/or SMASH_BUILD_DIR environment variables\n"
                        " are ignored when using SMASH_INSTALL_DIR to locate SMASH.\n")
    endif()
    if(NOT IS_ABSOLUTE $ENV{SMASH_INSTALL_DIR})
        message(FATAL_ERROR " \n" " SMASH_INSTALL_DIR must be set to a global path.\n")
    endif()
endif()

find_package(GSL 2.0 REQUIRED)
find_package(Eigen3 3.0 REQUIRED)
find_package(Pythia 8.316 EXACT REQUIRED)
set(SMASH_INCLUDE_DIR ${GSL_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIR} ${Pythia_INCLUDE_DIRS})
set(SMASH_LIBRARIES ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} ${Pythia_LIBRARIES} ${CMAKE_DL_LIBS})

if(DEFINED ENV{SMASH_INSTALL_DIR})

    if(DEFINED ENV{SMASH_VERSION})
        set(smash_folder "smash-$ENV{SMASH_VERSION}")
        if(NOT IS_DIRECTORY "$ENV{SMASH_INSTALL_DIR}/include/${smash_folder}")
            message(FATAL_ERROR " \n"
                                " No (complete) SMASH installation for version $ENV{SMASH_VERSION} found at\n"
                                " $ENV{SMASH_INSTALL_DIR}\n")
        endif()
    else()
        # Use 'include' installation folder to locate all installed versions, then consider last one
        file(GLOB smash_folders
             LIST_DIRECTORIES true
             RELATIVE "$ENV{SMASH_INSTALL_DIR}/include"
             "$ENV{SMASH_INSTALL_DIR}/include/smash-*")
        foreach(folder ${smash_folders})
            if(NOT IS_DIRECTORY "$ENV{SMASH_INSTALL_DIR}/include/${folder}")
                list(REMOVE_ITEM smash_folders "${folder}")
            endif()
        endforeach()
        list(LENGTH smash_folders number_of_versions)
        if(number_of_versions EQUAL 0)
            message(FATAL_ERROR " \n" " No SMASH installation found at $ENV{SMASH_INSTALL_DIR}\n")
        elseif(number_of_versions GREATER 1)
            message(STATUS "Multiple SMASH versions installed, CMake 3.18 needed to find last one")
            cmake_minimum_required(VERSION 3.18) # Required to sort versions in a natural way
            list(SORT smash_folders COMPARE NATURAL ORDER DESCENDING)
        endif()
        list(GET smash_folders 0 smash_folder)
        unset(number_of_versions)
        unset(smash_folders)
    endif()
    if(NOT IS_DIRECTORY "$ENV{SMASH_INSTALL_DIR}/lib/${smash_folder}")
        message(FATAL_ERROR " \n"
                            " Inconsistent installation of SMASH detected! Folder '${smash_folder}'\n"
                            " existing in 'include' but not found in 'lib' sub-folder\n"
                            " at $ENV{SMASH_INSTALL_DIR}\n")
    endif()
    # Now we can look for header files and libraries!
    list(APPEND
         SMASH_INCLUDE_DIR
         "$ENV{SMASH_INSTALL_DIR}/include/${smash_folder}/cuba"
         "$ENV{SMASH_INSTALL_DIR}/include/${smash_folder}/einhard"
         "$ENV{SMASH_INSTALL_DIR}/include/${smash_folder}/yaml-cpp"
         "$ENV{SMASH_INSTALL_DIR}/include/${smash_folder}")
    set(SMASH_INPUT_FILES_DIR "$ENV{SMASH_INSTALL_DIR}/share/${smash_folder}/input_files")
    set(path_for_libraries "$ENV{SMASH_INSTALL_DIR}/lib/${smash_folder}")
    find_library(SMASH_LIBRARY NAMES smash PATHS "${path_for_libraries}")
    find_library(EINHARD_LIBRARY NAMES einhard PATHS "${path_for_libraries}")
    find_library(CPPYAML_LIBRARY NAMES yaml-cpp PATHS "${path_for_libraries}")
    find_library(INTEGRATION_LIBRARY NAMES cuhre PATHS "${path_for_libraries}")
    unset(path_for_libraries)
    string(REPLACE "smash-" "" smash_found_version "${smash_folder}")

else() # Use SMASH_DIR and SMASH_BUILD_DIR --> DEPRECATED!

    list(APPEND
         SMASH_INCLUDE_DIR
         $ENV{SMASH_DIR}/3rdparty/Cuba-4.2.2
         $ENV{SMASH_DIR}/3rdparty/einhard
         $ENV{SMASH_DIR}/3rdparty/yaml-cpp-0.8.0/include
         $ENV{SMASH_DIR}/src/include
         ${SMASH_BUILD_DIR}/src/include # For the decaymodes and particles header files
    )
    set(SMASH_INPUT_FILES_DIR "$ENV{SMASH_DIR}/input")
    find_library(SMASH_LIBRARY NAMES smash PATHS ${SMASH_BUILD_DIR}/src)
    find_library(EINHARD_LIBRARY NAMES einhard PATHS ${SMASH_BUILD_DIR}/3rdparty/einhard)
    find_library(CPPYAML_LIBRARY NAMES yaml-cpp PATHS ${SMASH_BUILD_DIR}/3rdparty/yaml-cpp-0.8.0)
    find_library(INTEGRATION_LIBRARY NAMES cuhre
                 PATHS ${SMASH_BUILD_DIR}/3rdparty/Cuba-4.2.2/src/cuhre)
    set(smash_found_version "NOT-FOUND")

endif()

list(APPEND
     SMASH_LIBRARIES
     ${EINHARD_LIBRARY}
     ${CPPYAML_LIBRARY}
     ${INTEGRATION_LIBRARY}
     ${SMASH_LIBRARY})

# Validate selected include directories
foreach(folder ${SMASH_INCLUDE_DIR})
    if(NOT IS_DIRECTORY "${folder}")
        list(APPEND missing_folders "${folder}")
    endif()
endforeach()
list(LENGTH missing_folders number_of_missing_folders)
if(number_of_missing_folders GREATER 0)
    list(JOIN missing_folders "\n  - " missing_folders)
    message(FATAL_ERROR " \n" " Incomplete installation of SMASH. Missing folder(s):\n"
                        "  - ${missing_folders}\n")
endif()
unset(missing_folders)
unset(number_of_missing_folders)

# Validate found installation in standard CMake way
find_package_handle_standard_args(SMASH
                                  REQUIRED_VARS SMASH_LIBRARY
                                                EINHARD_LIBRARY
                                                CPPYAML_LIBRARY
                                                INTEGRATION_LIBRARY
                                                SMASH_INCLUDE_DIR
                                                SMASH_INPUT_FILES_DIR
                                  VERSION_VAR smash_found_version)
mark_as_advanced(SMASH_LIBRARY EINHARD_LIBRARY CPPYAML_LIBRARY INTEGRATION_LIBRARY)

if(SMASH_FOUND)
    # Some output displayed when --log-level=VERBOSE option is given to cmake command
    list(JOIN SMASH_INCLUDE_DIR "\n    " tmp_string)
    message(VERBOSE "SMASH includes:\n    ${tmp_string}")
    list(JOIN SMASH_LIBRARIES "\n    " tmp_string)
    message(VERBOSE "SMASH libraries:\n    ${tmp_string}")
    unset(tmp_string)
endif()

# Always inform user about outcome
message(STATUS "SMASH found: ${SMASH_FOUND}")
