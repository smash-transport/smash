# - Try to find ROOT
# This will define:
#
#  ROOT_FOUND
#  ROOT_CONFIG
#  ROOT_INCLUDE_DIR
#  ROOT_LIBRARIES
#  ROOT_VERSION

########################################################
#
#    Copyright (c) 2014-2015
#      SMASH Team
#
#    BSD 3-clause license
# 
######################################################### 

include(FindPackageHandleStandardArgs)
find_program(ROOT_CONFIG NAMES root-config HINTS $ENV{ROOTSYS}/bin)
if(ROOT_CONFIG)
   execute_process(COMMAND "${ROOT_CONFIG}" --prefix
      OUTPUT_VARIABLE ROOTSYS OUTPUT_STRIP_TRAILING_WHITESPACE)

   execute_process(COMMAND "${ROOT_CONFIG}" --version
      OUTPUT_VARIABLE ROOT_VERSION OUTPUT_STRIP_TRAILING_WHITESPACE)
   string(REPLACE "/" "." ROOT_VERSION "${ROOT_VERSION}")

   execute_process(COMMAND "${ROOT_CONFIG}" --incdir
      OUTPUT_VARIABLE ROOT_INCDIR OUTPUT_STRIP_TRAILING_WHITESPACE)
   set(ROOT_INCLUDE_DIR ${ROOT_INCDIR} CACHE PATH "path to ROOT includes")

   execute_process(COMMAND "${ROOT_CONFIG}" --glibs
      OUTPUT_VARIABLE ROOT_LINK_FLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
   set(ROOT_LIBRARIES ${ROOT_LINK_FLAGS} CACHE PATH "linker flags for ROOT")
endif()

FIND_PACKAGE_HANDLE_STANDARD_ARGS(ROOT
   REQUIRED_VARS ROOT_CONFIG ROOT_INCLUDE_DIR ROOT_LIBRARIES
   VERSION_VAR ROOT_VERSION
   FAIL_MESSAGE "required ROOT library not found"
   )

mark_as_advanced(ROOT_CONFIG)
