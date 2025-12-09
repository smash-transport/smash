########################################################
#
#    Copyright (c) 2025
#      SMASH Team
#
#    BSD 3-clause license
#
########################################################

# cmake-format: off
#=============================================================================
# This file was taken from the original Eigen3 repository and it has later
# been modified to support new library versions (when version 5.x has come out)
#=============================================================================
# - Try to find Eigen3 lib
#
# This module supports requiring a minimum version, e.g. you can do
#   find_package(Eigen3 3.1.2)
# to require version 3.1.2 or newer of Eigen3.
#
# Once done this will define
#
#  EIGEN3_FOUND - system has eigen lib with correct version
#  EIGEN3_INCLUDE_DIR - the eigen include directory
#  EIGEN3_VERSION - eigen version
#
# This module reads hints about search locations from
# the following environment variables:
#
# EIGEN3_ROOT
# EIGEN3_ROOT_DIR
#
# Copyright (c) 2006, 2007 Montel Laurent, <montel@kde.org>
# Copyright (c) 2008, 2009 Gael Guennebaud, <g.gael@free.fr>
# Copyright (c) 2009 Benoit Jacob <jacob.benoit.1@gmail.com>
# Redistribution and use is allowed according to the terms of the 2-clause BSD license.
#=============================================================================
# cmake-format: on

# ------------------ NEW: First try config mode (Eigen >= 5) -------------------

# Note that it is safe to use find_package here, and no infinite recursion occurs, because we are
# calling find_package() in CONFIG mode, which explicitly bypasses the module search. Furthermore,
# it is important to stress that the requested version here should be 5 and not 5.0 as the latter
# would reject versions 5.1, 5.2 and so on and this is not what we want to achieve. For the record,
# this module will need to be adjusted when Eigen3 version 6 or higher comes out. If the version
# requirement is dropped and simply find_package(Eigen3 QUIET CONFIG) is issued, then this will
# inherit the requirement from the parent scope and this is not wished. It can be patched unsetting
# and later restoring the Eigen3_FIND_VERSION_COMPLETE variable, but this has not yet been done
# because when version 6 of Eigen3 comes out it will probably possible to drop this entire file and
# simply rely on configuration mode finding requiring version 5 as minimum one.
find_package(Eigen3 5 QUIET CONFIG)

if(Eigen3_FOUND)
    if(NOT DEFINED Eigen3_VERSION AND NOT DEFINED Eigen3_VERSION_STRING)
        message(FATAL_ERROR " \n" " Eigen3 was found via CONFIG mode, but no version"
                            " information was provided.\n A valid Eigen3Config.cmake"
                            " must define Eigen3_VERSION or Eigen3_VERSION_STRING.\n")
    else()
        # Weird corner case, but it does not harm and allows later to just use Eigen3_VERSION
        if(NOT DEFINED Eigen3_VERSION)
            set(Eigen3_VERSION "${Eigen3_VERSION_STRING}")
        endif()
    endif()
    if(Eigen3_VERSION VERSION_LESS Eigen3_FIND_VERSION)
        message(FATAL_ERROR " \n " "Found Eigen3 version ${Eigen3_VERSION}, but "
                            "version ${Eigen3_FIND_VERSION} or newer is required.\n")
    else()
        # Map modern variables to legacy Find-module variables
        get_target_property(_eigen_inc Eigen3::Eigen INTERFACE_INCLUDE_DIRECTORIES)

        set(EIGEN3_INCLUDE_DIR "${_eigen_inc}")
        set(EIGEN3_VERSION "${Eigen3_VERSION}")
        set(EIGEN3_FOUND TRUE)
        mark_as_advanced(EIGEN3_INCLUDE_DIR)

        # Done — skip legacy detection
        message(STATUS "Eigen3 found via CONFIG mode (version ${EIGEN3_VERSION})")
        return()
    endif()
endif()

# -------------------- FALLBACK: Legacy Eigen 3.x detection --------------------

# Set the requested minimum to a proper version variable if not already the case
if(NOT Eigen3_FIND_VERSION)
    if(NOT Eigen3_FIND_VERSION_MAJOR)
        set(Eigen3_FIND_VERSION_MAJOR 2)
    endif(NOT Eigen3_FIND_VERSION_MAJOR)
    if(NOT Eigen3_FIND_VERSION_MINOR)
        set(Eigen3_FIND_VERSION_MINOR 91)
    endif(NOT Eigen3_FIND_VERSION_MINOR)
    if(NOT Eigen3_FIND_VERSION_PATCH)
        set(Eigen3_FIND_VERSION_PATCH 0)
    endif(NOT Eigen3_FIND_VERSION_PATCH)

    set(Eigen3_FIND_VERSION
        "${Eigen3_FIND_VERSION_MAJOR}.${Eigen3_FIND_VERSION_MINOR}.${Eigen3_FIND_VERSION_PATCH}")
endif(NOT Eigen3_FIND_VERSION)

macro(_eigen3_check_version)
    file(READ "${EIGEN3_INCLUDE_DIR}/Eigen/src/Core/util/Macros.h" _eigen3_version_header)

    string(REGEX MATCH "define[ \t]+EIGEN_WORLD_VERSION[ \t]+([0-9]+)" _eigen3_world_version_match
                 "${_eigen3_version_header}")
    set(EIGEN3_WORLD_VERSION "${CMAKE_MATCH_1}")
    string(REGEX MATCH "define[ \t]+EIGEN_MAJOR_VERSION[ \t]+([0-9]+)" _eigen3_major_version_match
                 "${_eigen3_version_header}")
    set(EIGEN3_MAJOR_VERSION "${CMAKE_MATCH_1}")
    string(REGEX MATCH "define[ \t]+EIGEN_MINOR_VERSION[ \t]+([0-9]+)" _eigen3_minor_version_match
                 "${_eigen3_version_header}")
    set(EIGEN3_MINOR_VERSION "${CMAKE_MATCH_1}")

    set(EIGEN3_VERSION ${EIGEN3_WORLD_VERSION}.${EIGEN3_MAJOR_VERSION}.${EIGEN3_MINOR_VERSION})
    if(${EIGEN3_VERSION} VERSION_LESS ${Eigen3_FIND_VERSION})
        set(EIGEN3_VERSION_OK FALSE)
    else(${EIGEN3_VERSION} VERSION_LESS ${Eigen3_FIND_VERSION})
        set(EIGEN3_VERSION_OK TRUE)
    endif(${EIGEN3_VERSION} VERSION_LESS ${Eigen3_FIND_VERSION})

    if(NOT EIGEN3_VERSION_OK)

        message(STATUS "Eigen3 version ${EIGEN3_VERSION} found in ${EIGEN3_INCLUDE_DIR}, "
                       "but at least version ${Eigen3_FIND_VERSION} is required")
    endif(NOT EIGEN3_VERSION_OK)
endmacro(_eigen3_check_version)

if(EIGEN3_INCLUDE_DIR)

    # in cache already
    _eigen3_check_version()
    set(EIGEN3_FOUND ${EIGEN3_VERSION_OK})

    # Abort if the version could not be extracted or wasn't matching requirements
    if(NOT EIGEN3_FOUND)
        message(FATAL_ERROR "Unsuitable Eigen3 version '${EIGEN3_VERSION}' found.")
    endif()

else(EIGEN3_INCLUDE_DIR)

    find_path(EIGEN3_INCLUDE_DIR
              NAMES signature_of_eigen3_matrix_library
              HINTS ENV EIGEN3_ROOT ENV EIGEN3_ROOT_DIR
              PATHS ${CMAKE_INSTALL_PREFIX}/include ${KDE4_INCLUDE_DIR}
              PATH_SUFFIXES eigen3 eigen)

    if(EIGEN3_INCLUDE_DIR)
        _eigen3_check_version()
    endif(EIGEN3_INCLUDE_DIR)

    include(FindPackageHandleStandardArgs)
    find_package_handle_standard_args(Eigen3 DEFAULT_MSG EIGEN3_INCLUDE_DIR EIGEN3_VERSION_OK)

    mark_as_advanced(EIGEN3_INCLUDE_DIR)

endif(EIGEN3_INCLUDE_DIR)
