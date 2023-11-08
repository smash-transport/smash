#!/usr/bin/env bash

if [[ $# -ne 1 ]] && [[ $# -ne 2 ]]; then
  echo "Script requires 2 or 3 arguments. Usage: $0 SMASH_DIR [PYTHIA_CONFIG_PATH]"
  echo "PYTHIA_CONFIG_PATH needs to be given if Pythia is not installed system wide."
  exit 1
fi

SMASH_DIR_VAR=$1
PYTHIA_CONFIG_PATH=$2

if [[ "${PYTHIA_CONFIG_PATH}" = '' ]]; then
    PYTHIA_OPTION_FOR_CMAKE=''
else
    PYTHIA_OPTION_FOR_CMAKE="-DPythia_CONFIG_EXECUTABLE=${PYTHIA_CONFIG_PATH}"
fi

export SMASH_DIR="${SMASH_DIR_VAR}"
EXAMPLE_DIR=${SMASH_DIR_VAR}/examples/using_SMASH_as_library

fail_and_rm_build() {
  echo "$*" >&2
  cd .. && rm -r build
  exit 1
}

cd "${EXAMPLE_DIR}"
mkdir -p build && cd build
# Assume the default compiler is fine or that `CC` and `CXX` environment variables are set
cmake "${PYTHIA_OPTION_FOR_CMAKE}" .. || fail_and_rm_build "Failed to setup cmake for SMASH library examples"
make -j$(nproc) || fail_and_rm_build "Failed to build SMASH library examples"

./example || fail_and_rm_build "Failed to run SMASH library examples"
./example_rate_equations || fail_and_rm_build "Failed to run SMASH rate equation library examples"
mkdir data  # create output directory for smash wrapper
./example_smash_wrapper || fail_and_rm_build "Failed to execute SMASH wrapper library example"

cd .. && rm -r build # Clean-up
