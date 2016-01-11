#############################################################
# Cmake code to integrate Pythia contributed by 
# K. Gallmeister, Goethe University, March 2015
# <gallmei@th.physik.uni-frankfurt.de>
###########################################################

# Try to find a Pythia 8.x installation
#
# Relies on info given by the executable 'pythia8-config'
#
# This package defines
#  Pythia_FOUND - Pythia has been found
#  Pythia_INCLUDE_DIRS - The directory in which the Pythia headers reside
#  Pythia_LIBRARIES - The Pythia library / libraries
#  Pythia_LHAPDFDummy_LIBRARY - The LHAPDF dummy library, use if the real LHAPDF is not installed / found


FIND_PROGRAM(Pythia_CONFIG_EXECUTABLE NAMES pythia8-config HINTS ../3rdparty/pythia8215/bin)
IF(${Pythia_CONFIG_EXECUTABLE} MATCHES "Pythia_CONFIG_EXECUTABLE-NOTFOUND")
  MESSAGE(STATUS "Looking for Pythia... - pythia8-config executable not found")
ELSE(${Pythia_CONFIG_EXECUTABLE} MATCHES "Pythia_CONFIG_EXECUTABLE-NOTFOUND")
  MESSAGE(STATUS "Looking for Pythia... - using pythia8-config executable")
  EXEC_PROGRAM(${Pythia_CONFIG_EXECUTABLE} ARGS "--libdir" OUTPUT_VARIABLE Pythia_LIBDIR)
  FIND_LIBRARY( Pythia_LIBRARY
    NAMES pythia8
    PATHS ${Pythia_LIBDIR}
    PATH_SUFFIXES PYTHIA pythia Pythia PYTHIA8 Pythia8 pythia8  # suggest some path suffixes in which the headers could be located
    )
 # FIND_LIBRARY( Pythia_LHAPDFDummy_LIBRARY
 #   NAMES lhapdfdummy
 #   PATHS ${Pythia_LIBDIR}
 #   PATH_SUFFIXES PYTHIA pythia Pythia PYTHIA8 Pythia8 pythia8  # suggest some path suffixes in which the headers could be located
 # )
  EXEC_PROGRAM(${Pythia_CONFIG_EXECUTABLE} ARGS "--includedir" OUTPUT_VARIABLE Pythia_INCLUDE_DIR)
  EXEC_PROGRAM(${Pythia_CONFIG_EXECUTABLE} ARGS "--xmldoc" OUTPUT_VARIABLE Pythia_xmldoc_PATH)
ENDIF()


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
ENDIF( Pythia_FOUND )

# the variables will only show up in the GUI in the "advanced" view
MARK_AS_ADVANCED(Pythia_INCLUDE_DIR Pythia_LIBRARY Pythia_LHAPDFDummy_LIBRARY Pythia_xmldoc_PATH)
MARK_AS_ADVANCED(Pythia_CONFIG_EXECUTABLE)
