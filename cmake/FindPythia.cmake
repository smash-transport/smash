#############################################################
# Cmake code to integrate Pythia contributed by 
# K. Gallmeister, Goethe University, March 2015
# <gallmei@th.physik.uni-frankfurt.de>
###########################################################

# Try to find a Pythia 8.x installation
#
# Will use info given by the executable 'pythia8-config'. Otherwise, tries to use a local build.
#
# This package defines
#  Pythia_FOUND - Pythia has been found
#  Pythia_INCLUDE_DIRS - The directory in which the Pythia headers reside
#  Pythia_LIBRARIES - The Pythia library / libraries
#  Pythia_LHAPDFDummy_LIBRARY - The LHAPDF dummy library, use if the real LHAPDF is not installed / found


FIND_PROGRAM(Pythia_CONFIG_EXECUTABLE NAMES pythia8-config)
IF(${Pythia_CONFIG_EXECUTABLE} MATCHES "Pythia_CONFIG_EXECUTABLE-NOTFOUND")
  MESSAGE(STATUS "Looking for Pythia... - pythia8-config executable not found, please specify it with \"-DPythia_CONFIG_EXECUTABLE\"")
  #SET( Pythia_INCLUDE_DIR ../3rdparty/pythia8230/include )
  #SET( Pythia_LIBDIR ../3rdparty/pythia8230/lib )
  #SET( Pythia_xmldoc_PATH ../3rdparty/pythia8230/share/Pythia8/xmldoc )
ELSE()
  MESSAGE(STATUS "Looking for Pythia... - using pythia8-config executable")
  EXEC_PROGRAM(${Pythia_CONFIG_EXECUTABLE} ARGS "--libdir" OUTPUT_VARIABLE Pythia_LIBDIR)
  EXEC_PROGRAM(${Pythia_CONFIG_EXECUTABLE} ARGS "--includedir" OUTPUT_VARIABLE Pythia_INCLUDE_DIR)
  EXEC_PROGRAM(${Pythia_CONFIG_EXECUTABLE} ARGS "--xmldoc" OUTPUT_VARIABLE Pythia_xmldoc_PATH)
  EXEC_PROGRAM(${Pythia_CONFIG_EXECUTABLE} ARGS "--prefix" OUTPUT_VARIABLE Pythia_prefix)
  FILE(READ ${Pythia_INCLUDE_DIR}/Pythia8/Pythia.h Pythia_header)
  STRING(REGEX MATCH "#define PYTHIA_VERSION_INTEGER ([0-9])([0-9][0-9][0-9])" _ ${Pythia_header})
  SET(Pythia_VERSION_MAJOR ${CMAKE_MATCH_1})
  SET(Pythia_VERSION_MINOR ${CMAKE_MATCH_2})
  SET(Pythia_VERSION ${Pythia_VERSION_MAJOR}.${Pythia_VERSION_MINOR})
  if(NOT(${Pythia_VERSION} VERSION_EQUAL ${Pythia_FIND_VERSION}))
    MESSAGE( STATUS "** WRONG Pythia version: ${Pythia_VERSION},"
                    " required ${Pythia_FIND_VERSION}." )
    set(Pythia_VERSION_OK FALSE)
  else(${Pythia_VERSION} VERSION_LESS ${Pythia_FIND_VERSION})
    MESSAGE( STATUS "** Pythia version ok: ${Pythia_VERSION},"
                    " required ${Pythia_FIND_VERSION}.")
    set(Pythia_VERSION_OK TRUE)
  endif()
ENDIF()

FIND_LIBRARY( Pythia_LIBRARY
  NAMES pythia8
  PATHS ${Pythia_LIBDIR}
  PATH_SUFFIXES PYTHIA pythia Pythia PYTHIA8 Pythia8 pythia8  # suggest some path suffixes in which the headers could be located
)

# Check that the actual found library is in the same directory as the config executable
# first get directory from the lib-path 
get_filename_component(Pythia_actual_dir ${Pythia_LIBRARY} DIRECTORY)

# strip lib-component from actual found directory path
string(FIND ${Pythia_actual_dir} "lib" Pythia_actual_dir_lib_pos REVERSE)

string(SUBSTRING ${Pythia_actual_dir} 0 ${Pythia_actual_dir_lib_pos} Pythia_actual_lib_abs_path)

# Now the same for the used config
get_filename_component(Pythia_CONFIG_DIR ${Pythia_CONFIG_EXECUTABLE} DIRECTORY)
string(FIND ${Pythia_CONFIG_DIR} "bin" Pythia_CONFIG_DIR_bin_pos REVERSE)
string(SUBSTRING ${Pythia_CONFIG_DIR} 0 ${Pythia_CONFIG_DIR_bin_pos} Pythia_CONFIG_DIR_abs_path)

if(NOT (${Pythia_CONFIG_DIR_abs_path} STREQUAL ${Pythia_actual_lib_abs_path}))
  MESSAGE(WARNING "Path to Pythia library and config do not match. Consider cleaning your build-folder and run CMAKE again.")
endif()

# FIND_LIBRARY( Pythia_LHAPDFDummy_LIBRARY
#   NAMES lhapdfdummy
#   PATHS ${Pythia_LIBDIR}
#   PATH_SUFFIXES PYTHIA pythia Pythia PYTHIA8 Pythia8 pythia8  # suggest some path suffixes in which the headers could be located
#   HINTS ../3rdparty/pythia8125/lib
# )


# adhere to the standard nomenclature for find_package routines and set some variables
SET( Pythia_INCLUDE_DIRS ${Pythia_INCLUDE_DIR} )
SET( Pythia_LIBRARIES ${Pythia_LIBRARY} )

# handle the QUIETLY and REQUIRED arguments and set Pythia_FOUND to TRUE
# if all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
#FIND_PACKAGE_HANDLE_STANDARD_ARGS( Pythia DEFAULT_MSG Pythia_LIBRARY Pythia_INCLUDE_DIR Pythia_xmldoc_PATH Pythia_LHAPDFDummy_LIBRARY)
FIND_PACKAGE_HANDLE_STANDARD_ARGS( Pythia DEFAULT_MSG Pythia_LIBRARY Pythia_INCLUDE_DIR Pythia_xmldoc_PATH)
SET( Pythia_FOUND ${PYTHIA_FOUND} CACHE INTERNAL "Provide Pythia_FOUND in addition to PYTHIA_FOUND" FORCE )


# display some status information
IF( Pythia_FOUND )
  MESSAGE( STATUS "** Pythia 8 library: ${Pythia_LIBRARIES}" )
  MESSAGE( STATUS "** Pythia 8 include: ${Pythia_INCLUDE_DIRS}" )
  MESSAGE( STATUS "** Pythia 8 xmldoc:  ${Pythia_xmldoc_PATH}" )
ENDIF()
IF (NOT Pythia_VERSION_OK)
  SET( Pythia_FOUND FALSE)
ENDIF ()

# the variables will only show up in the GUI in the "advanced" view
MARK_AS_ADVANCED(Pythia_INCLUDE_DIR Pythia_LIBRARY Pythia_LHAPDFDummy_LIBRARY Pythia_xmldoc_PATH)
MARK_AS_ADVANCED(Pythia_CONFIG_EXECUTABLE)
