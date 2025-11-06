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
# This will define:
#
#  RIVET_CONFIG
#  RIVET_VERSION
#  RIVET_INCLUDE_DIR
#  RIVET_LIBRARIES
#  Rivet_FOUND
#=======================================================
# cmake-format: on

find_program(RIVET_CONFIG rivet-config DOC "Path to the Rivet configuration script")

if(RIVET_CONFIG)
    execute_process(COMMAND ${RIVET_CONFIG} --version OUTPUT_VARIABLE RIVET_VERSION
                    OUTPUT_STRIP_TRAILING_WHITESPACE)
    execute_process(COMMAND ${RIVET_CONFIG} --includedir OUTPUT_VARIABLE RIVET_INCLUDE_DIR
                    OUTPUT_STRIP_TRAILING_WHITESPACE)
    execute_process(COMMAND ${RIVET_CONFIG} --libs OUTPUT_VARIABLE RIVET_LIBRARIES
                    OUTPUT_STRIP_TRAILING_WHITESPACE)
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Rivet REQUIRED_VARS RIVET_CONFIG RIVET_INCLUDE_DIR RIVET_LIBRARIES
                                  VERSION_VAR RIVET_VERSION)

mark_as_advanced(RIVET_CONFIG)
