########################################################
#
#    Copyright (c) 2015-2016,2018-2022,2024-2025
#      SMASH Team
#
#    BSD 3-clause license
#
########################################################

# cmake-format: off
#=============================================================================
# This cmake code to integrate Pythia contains contributions by
#   K. Gallmeister, Goethe University, March 2015
#   <gallmei@th.physik.uni-frankfurt.de>
# and by
#   A. Verbytskyi, Max-Planck Institute für Physics, January 2022
#   <andrii.verbytskyi@mpp.mpg.de>
#=============================================================================
# - Locate pythia library
#
# Exploits the environment variables:
#  PYTHIA_ROOT_DIR or PYTHIA8
# or the options
#  -DPythia_CONFIG_EXECUTABLE or -DPYTHIA_ROOT_DIR
#
# Tries to find the version defined in the variable:
#
#  Pythia_VERSION
#
# Defines:
#
#  Pythia_FOUND
#  Pythia_VERSION
#  Pythia_INCLUDE_DIR
#  PYTHIA_XMLDOC_DIR
#  Pythia_xmldoc_PATH
#  Pythia_INCLUDE_DIRS
#  Pythia_LIBRARY
#  Pythia_LIBDIR
#  Pythia_LHAPDFDummy_LIBRARY
#  Pythia_LIBRARIES : includes 3 libraries above; not to be used if lhapdf is used
#=============================================================================
# cmake-format: on

if(NOT ("${Pythia_CONFIG_EXECUTABLE}" STREQUAL ""))
    message(STATUS "Trying to locate Pythia using variable Pythia_CONFIG_EXECUTABLE = ${Pythia_CONFIG_EXECUTABLE}"
    )
    find_program(Pythia_CONFIG_EXECUTABLE NAMES pythia8-config)
    if(${Pythia_CONFIG_EXECUTABLE} MATCHES "Pythia_CONFIG_EXECUTABLE-NOTFOUND")
        message(FATAL_ERROR "pythia8-config executable not found, please check \"-DPythia_CONFIG_EXECUTABLE\" or use \"-DPYTHIA_ROOT_DIR\""
        )
    else()
        execute_process(COMMAND ${Pythia_CONFIG_EXECUTABLE} --prefix OUTPUT_VARIABLE PYTHIA_ROOT_DIR
                        OUTPUT_STRIP_TRAILING_WHITESPACE)
    endif()
endif()

if("${PYTHIA_ROOT_DIR}" STREQUAL "")
    if(DEFINED ENV{PYTHIA8})
        message(STATUS "Trying to locate Pythia using the environment variable PYTHIA8 = $ENV{PYTHIA8}"
        )
        set(PYTHIA_ROOT_DIR $ENV{PYTHIA8})
    elseif(DEFINED ENV{PYTHIA_ROOT_DIR})
        message(STATUS "Trying to locate Pythia using the environment variable PYTHIA_ROOT_DIR = $ENV{PYTHIA_ROOT_DIR}"
        )
        set(PYTHIA_ROOT_DIR $ENV{PYTHIA_ROOT_DIR})
    else()
        message(STATUS "The installation directory of Pythia can be specified by setting the environment variables PYTHIA_ROOT_DIR or PYTHIA8\n"
                       "   or either with \"-DPythia_CONFIG_EXECUTABLE\" or \"-DPYTHIA_ROOT_DIR\", but none of these options was used.\n"
                       "   Trying to proceed assuming that Pythia is installed under /usr.")
        set(PYTHIA_ROOT_DIR "/usr")
    endif()
endif()

find_path(Pythia_INCLUDE_DIR Pythia.h Pythia8/Pythia.h HINTS ${PYTHIA_ROOT_DIR}/include)

find_path(Pythia_XMLDOC_DIR Version.xml
          HINTS ${PYTHIA_ROOT_DIR}/xmldoc ${PYTHIA_ROOT_DIR}/share/Pythia8/xmldoc
                ${PYTHIA_ROOT_DIR}/share/pythia8-data/xmldoc
                ${PYTHIA_ROOT_DIR}/share/doc/packages/pythia/xmldoc)

if(Pythia_INCLUDE_DIR AND Pythia_XMLDOC_DIR)
    find_library(Pythia_LIBRARY NAMES pythia8 Pythia8 HINTS ${PYTHIA_ROOT_DIR}/lib
                                                            ${PYTHIA_ROOT_DIR}/lib64)
    get_filename_component(Pythia_LIBDIR ${Pythia_LIBRARY} DIRECTORY)
    find_library(Pythia_LHAPDFDummy_LIBRARY NAMES lhapdfdummy HINTS ${PYTHIA_ROOT_DIR}/lib
                                                                    ${PYTHIA_ROOT_DIR}/lib64)
    set(Pythia_INCLUDE_DIRS ${Pythia_INCLUDE_DIR} ${Pythia_INCLUDE_DIR}/Pythia
                            ${Pythia_INCLUDE_DIR}/PythiaPlugins)
    set(Pythia_LIBRARIES ${Pythia_LIBRARY})
    file(READ ${Pythia_INCLUDE_DIR}/Pythia8/Pythia.h Pythia_header)
    string(REGEX MATCH "#define PYTHIA_VERSION_INTEGER ([0-9])([0-9][0-9][0-9])" _ ${Pythia_header})
    set(Pythia_VERSION_MAJOR ${CMAKE_MATCH_1})
    set(Pythia_VERSION_MINOR ${CMAKE_MATCH_2})
    set(Pythia_VERSION ${Pythia_VERSION_MAJOR}.${Pythia_VERSION_MINOR})
else()
    set(Pythia_FOUND FALSE)
endif()

set(Pythia_INCLUDE_DIRS ${Pythia_INCLUDE_DIR})
set(Pythia_xmldoc_PATH ${Pythia_XMLDOC_DIR})

# handle the QUIETLY and REQUIRED arguments
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Pythia REQUIRED_VARS Pythia_LIBRARY Pythia_INCLUDE_DIRS
                                                       Pythia_xmldoc_PATH
                                  VERSION_VAR Pythia_VERSION)

# display some status information
if(Pythia_FOUND)
    message(STATUS "Pythia 8 library: ${Pythia_LIBRARIES}")
    message(STATUS "Pythia 8 include: ${Pythia_INCLUDE_DIRS}")
    message(STATUS "Pythia 8 xmldoc:  ${Pythia_xmldoc_PATH}")
endif()

# the variables listed here will only show up in the GUI (ccmake) in the "advanced" view
mark_as_advanced(Pythia_FOUND
                 Pythia_INCLUDE_DIR
                 Pythia_LIBRARY
                 Pythia_LIBRARIES
                 Pythia_LHAPDFDummy_LIBRARY
                 Pythia_XMLDOC_DIR
                 Pythia_LIBDIR
                 Pythia_CONFIG_EXECUTABLE)
