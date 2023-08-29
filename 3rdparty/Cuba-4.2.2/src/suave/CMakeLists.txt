########################################################
#
#    Copyright (c) 2017,2021-2022
#      SMASH Team
#
#    BSD 3-clause license
#
#########################################################

include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/../common)

add_library(suave STATIC Suave.c ${COMMON_DEPS})
