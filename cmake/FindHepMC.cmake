# - Try to find HepMC
# This will define:
#
#  HEPMC_FOUND
#  HEPMC_INCLUDE_DIR
#  HEPMC_LIBRARIES

########################################################
#
#    Copyright (c) 2020
#      SMASH Team
#
#    BSD 3-clause license
#
#########################################################


find_path(HEPMC_INCLUDE_DIR NAMES HepMC3/Version.h
  PATHS
  ${HEPMC_DIR}/include/
  /usr/local/lib/include/
  /usr/local/include/
)

find_library(HEPMC_LIBRARIES NAMES libHepMC3_static.a libHepMC3.dylib libHepMC3.so
  PATHS
  ${HEPMC_DIR}/lib
  ${HEPMC_DIR}/lib64
  /usr/local/lib
  /usr/local/lib/lib64
)


include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(HepMC
  DEFAULT_MSG
  HEPMC_LIBRARIES
  HEPMC_INCLUDE_DIR)
