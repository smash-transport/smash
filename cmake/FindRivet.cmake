########################################################
#
#    Copyright (c) 2022,2025
#      SMASH Team
#
#    BSD 3-clause license
#
########################################################

# cmake-format: off
#=======================================================
# - Try to find Rivet
#
# This module defines:
#
#  RIVET_CONFIG          Path to rivet-config executable
#  RIVET_VERSION         Rivet version string
#  RIVET_INCLUDE_DIR     Rivet include directory
#  RIVET_LIBDIR          Rivet library directory
#  RIVET_LIBRARY         Rivet library file path
#  RIVET_LIBRARIES       Linker flags from rivet-config --libs
#  Rivet_FOUND           Whether Rivet was found
#
# It also defines the imported target Rivet::Rivet which provides:
#   - SYSTEM include directories
#   - INTERFACE link libraries and flags from rivet-config
#=======================================================
# cmake-format: on

find_program(RIVET_CONFIG rivet-config DOC "Path to the Rivet configuration script")

if(RIVET_CONFIG)
    execute_process(COMMAND ${RIVET_CONFIG} --version OUTPUT_VARIABLE RIVET_VERSION
                    OUTPUT_STRIP_TRAILING_WHITESPACE)
    execute_process(COMMAND ${RIVET_CONFIG} --includedir OUTPUT_VARIABLE RIVET_INCLUDE_DIR
                    OUTPUT_STRIP_TRAILING_WHITESPACE)
    execute_process(COMMAND ${RIVET_CONFIG} --libdir OUTPUT_VARIABLE RIVET_LIBDIR
                    OUTPUT_STRIP_TRAILING_WHITESPACE)
    execute_process(COMMAND ${RIVET_CONFIG} --libs OUTPUT_VARIABLE RIVET_LIBRARIES
                    OUTPUT_STRIP_TRAILING_WHITESPACE)

    # Convert whitespace-separated flags into a proper CMake list. This ensures consistency with how
    # CMake find-modules usually return lists of libraries and include directories.
    separate_arguments(RIVET_LIBRARIES NATIVE_COMMAND "${RIVET_LIBRARIES}")
endif()

if(RIVET_LIBDIR)
    find_library(RIVET_LIBRARY
                 NAMES Rivet
                 HINTS ${RIVET_LIBDIR}
                 NO_DEFAULT_PATH)
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Rivet REQUIRED_VARS RIVET_CONFIG RIVET_INCLUDE_DIR RIVET_LIBRARY
                                  VERSION_VAR RIVET_VERSION)

mark_as_advanced(RIVET_CONFIG RIVET_VERSION RIVET_LIBDIR RIVET_LIBRARIES)

if(Rivet_FOUND AND NOT TARGET Rivet::Rivet)
    add_library(Rivet::Rivet UNKNOWN IMPORTED)
    set_target_properties(Rivet::Rivet PROPERTIES IMPORTED_LOCATION "${RIVET_LIBRARY}")
    target_include_directories(Rivet::Rivet SYSTEM INTERFACE ${RIVET_INCLUDE_DIR})
endif()
