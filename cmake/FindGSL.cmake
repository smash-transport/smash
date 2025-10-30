########################################################
#
#    Copyright (c) 2014,2018,2022
#      SMASH Team
#
#    BSD 3-clause license
#
########################################################

# cmake-format: off
#=============================================================================
# - Try to find GSL (GNU Scientific Library)
#
# This will define:
#
#  GSL_FOUND
#  GSL_INCLUDE_DIR
#  GSL_LIBRARY
#  GSL_CBLAS_LIBRARY
#
#  This module is based on the module provided by cmake starting with version 3.2
#  See https://cmake.org/licensing for copyright info.
#
#=============================================================================
# cmake-format: on

include(FindPackageHandleStandardArgs)

# first check if GSL_ROOT_DIR is set (either as environment variable or supplied as cmake-option. If
# so, use it
if(EXISTS "$ENV{GSL_ROOT_DIR}")
    file(TO_CMAKE_PATH "$ENV{GSL_ROOT_DIR}" GSL_ROOT_DIR)
    set(GSL_ROOT_DIR "${GSL_ROOT_DIR}" CACHE PATH "Prefix for GSL installation")

elseif(EXISTS "${GSL_ROOT_DIR}")
    file(TO_CMAKE_PATH ${GSL_ROOT_DIR} GSL_ROOT_DIR)
    set(GSL_ROOT_DIR "${GSL_ROOT_DIR}" CACHE PATH "Prefix for GSL installation")
endif()

# no user supplied path. Try to find gsl on our own.
if(NOT EXISTS "${GSL_ROOT_DIR}")
    set(GSL_USE_PKGCONFIG ON)
endif()

if(GSL_USE_PKGCONFIG)
    find_package(PkgConfig)
    pkg_check_modules(GSL gsl)
    if(EXISTS "${GSL_INCLUDE_DIR}")
        get_filename_component(GSL_ROOT_DIR "${GSL_INCLUDE_DIR}" PATH CACHE)
    endif()
endif()

find_path(GSL_INCLUDE_DIR NAMES gsl HINTS ${GSL_ROOT_DIR}/include ${GSL_INCLUDE_DIR})

find_library(GSL_LIBRARY NAMES gsld gsl HINTS ${GSL_ROOT_DIR}/lib ${GSL_LIBDIR})

find_library(GSL_CBLAS_LIBRARY NAMES gslcblas cblas HINTS ${GSL_ROOT_DIR}/lib ${GSL_LIBDIR})

set(GSL_INCLUDE_DIRS ${GSL_INCLUDE_DIR})
set(GSL_LIBRARIES ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY})

# we could probably also search in the filepath for the version.
if(NOT GSL_VERSION)
    find_program(GSL_CONFIG_EXE NAMES gsl-config HINTS ${GSL_ROOT_DIR}/bin ${GSL_ROOT_DIR})
    if(EXISTS ${GSL_CONFIG_EXE})
        execute_process(COMMAND "${GSL_CONFIG_EXE}" "--version"
                        RESULT_VARIABLE GSL_PROC_STATUS
                        OUTPUT_VARIABLE GSL_VERSION
                        OUTPUT_STRIP_TRAILING_WHITESPACE)
    endif()
endif()

find_package_handle_standard_args(GSL
                                  FOUND_VAR GSL_FOUND
                                  REQUIRED_VARS GSL_INCLUDE_DIR GSL_LIBRARY GSL_CBLAS_LIBRARY
                                  VERSION_VAR GSL_VERSION)

mark_as_advanced(GSL_CBLAS_LIBRARY GSL_CONFIG_EXE GSL_INCLUDE_DIRS GSL_LIBRARIES)
