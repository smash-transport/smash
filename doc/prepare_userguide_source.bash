#!/usr/bin/env bash
#
########################################################
#
#    Copyright (c) 2014,2022
#      SMASH Team
#
#     GNU General Public License (GPLv3 or later)
#
########################################################
#
# This script, generally called by CMake when building the user guide,
# is supposed to collect all the /*!\Userguide documentation blocks into
# into a single file.

# Parse command line arguments. Enforce usage of '--' separator between
# command line options and arguments which refer to filenames.
if [[ $3 != '--' ]]; then
    printf "\e[91m-- Script '$(basename ${BASH_SOURCE[0]})' not properly called.\e[0m\n"
    exit 1
fi
readonly sed_executable="$1"
readonly output_filename="$2"
shift 3 # Ignore -- argument, too.

# Collect user guide documentation blocks. Each block is meant to be
# a Doxygen documentation block either starting by /** or by /*! and
# necessarily ending with */ (other Doxygen styles as described in the
# official guide https://doxygen.nl/manual/docblocks.html are not
# accepted). On the starting block line, the Doxygen alias \Userguide
# has to be found, optionally separated by spaces from the block initial
# delimiter. Further optional accepted, although discouraged, content:
#   - On the starting line of the block, characters after \Userguide
#     can appear, provided that there is some space after the alias.
#   - On the ending line of the block, space can appear before and/or
#     after the delimiter, but nothing else can occur on that line.
#
# Note that the use of -E is needed to use (a|b) regex, which would
# still be interpreted and work on GNU sed, but not on FreeBSD or
# macOS. For this particular case, so far, no distinction between
# the two is needed, which might be the case if a word boundary was
# needed (this is \< or \> in GNU sed and [[:<:]] or [[:>:]] on macOS).

${sed_executable} -E '/[/][*][*!][[:space:]]*[\]Userguide($|[[:space:]])/,/^[[:space:]]*[*][/][[:space:]]*$/!d' "$@" > "${output_filename}"
if [[ $? -ne 0 || ! -s "${output_filename}" ]]; then
    printf "\e[91m-- Error occurred collecting user guide blocks!\e[0m\n"
    exit 1
fi
