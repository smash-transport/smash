#!/usr/bin/env bash

########################################################
#
#    Copyright (c) 2015,2019-2020,2022,2026
#      SMASH Team
#
#    GNU General Public License (GPLv3 or later)
#
########################################################

# Stricter bash mode
set -euo pipefail

trap 'printf "\n"' EXIT

function fail()
{
    local -r label='ERROR: '
    printf "\n \e[1;91m${label}\e[22m$1\e[0m\n" 1>&2
    shift
    if [[ $# -gt 0 ]]; then
        printf " \e[91m${label//?/ }%s\e[0m\n" "$@" 1>&2
    fi
    exit 1
}

if [[ $# -ne 2 ]] && [[ $# -ne 3 ]]; then
  fail "Script requires 2 or 3 command line arguments." \
       "Usage: $0 SMASH_DIR SMASH_BUILD_DIR [PYTHIA_CONFIG_PATH]" \
       'PYTHIA_CONFIG_PATH needs to be given if Pythia is not installed system wide.'
fi

SMASH_DIR_VAR=$1
SMASH_BUILD_DIR=$2
PYTHIA_CONFIG_PATH=$3
EXAMPLE_DIR=${SMASH_DIR_VAR}/examples/using_SMASH_as_library

# Few sanity checks.
if [[ ! -d "${EXAMPLE_DIR}" ]]; then
  fail "Path to SMASH codebase invalid. Check the arguments passed to the script."
elif [[ ! -d "${SMASH_BUILD_DIR}" ]]; then
  fail "Path to SMASH build directory invalid. Check the arguments passed to the script."
fi

if [[ "${PYTHIA_CONFIG_PATH}" = '' ]]; then
    PYTHIA_OPTION_FOR_CMAKE=''
else
    PYTHIA_OPTION_FOR_CMAKE="-DPythia_CONFIG_EXECUTABLE=${PYTHIA_CONFIG_PATH}"
fi

function do_clean_up()
{
  printf 'Removing folders:\n  %s\n  %s\n' "${EXAMPLE_DIR}/build" "${SMASH_INSTALL_DIR}"
  rm -rf "${EXAMPLE_DIR}/build" "${SMASH_INSTALL_DIR}"
}

function fail_and_rm_build() {
  do_clean_up
  fail "$@"
}

# SMASH needs to be installed before using it as a library. We do this in a temporary folder.
# However, we need to run CMake again in the build folder to be sure to set the correct installation
# prefix. This should work both if SMASH was not set up and even if it was already built.
SMASH_INSTALL_DIR="$(mktemp -d)"
cd "${SMASH_BUILD_DIR}"
cmake -DCMAKE_INSTALL_PREFIX="${SMASH_INSTALL_DIR}"\
      "${PYTHIA_OPTION_FOR_CMAKE}" "${SMASH_DIR_VAR}" || fail_and_rm_build 'Failed to set up SMASH'
make -j$(nproc) install || fail_and_rm_build 'Failed to install SMASH'

# Now set up "external" CMake project and build it using SMASH as library
cd "${EXAMPLE_DIR}"
mkdir -p build && cd build
# Assume the default compiler is fine or that `CC` and `CXX` environment variables are set
cmake -DCMAKE_PREFIX_PATH="${SMASH_INSTALL_DIR}" .. || fail_and_rm_build 'Failed to setup cmake for SMASH library examples'
make -j$(nproc) || fail_and_rm_build 'Failed to build SMASH library examples'

./example || fail_and_rm_build 'Failed to run SMASH library example'
./example_rate_equations || fail_and_rm_build 'Failed to run SMASH rate equation library example'
./example_smash_wrapper || fail_and_rm_build 'Failed to execute SMASH wrapper library example'

do_clean_up
