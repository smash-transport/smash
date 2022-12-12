#!/bin/sh
#
########################################################
#
#    Copyright (c) 2014,2022
#      SMASH Team
#
#    BSD 3-clause license
#
#########################################################
#
# This script, generally called by CMake when building
# the user guide, is supposed to collect into a single
# file all the /*!\Userguide documentation blocks.
SED="$1"
output="$2"
shift 2
$SED -e '/\/\*[*!] *\\Userguide\>/,/\*\//!d' "$@" > $output
