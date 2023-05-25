#!/usr/bin/env bash
# run with: bash run_dyIC.bash [SMASH_PATH] [EIGEN3_PATH] [PYTHIA8309_PATH]/bin/pythia8-config 
if [[ $# -ne 2 ]] && [[ $# -ne 3 ]]; then
  echo "Script requires 2 or 3 arguments. Usage: $0 SMASH_DIR EIGEN_DIR [PYTHIA_DIR_VAR]"
  echo "PYTHIA_DIR_VAR needs to be given if Pythia is not installed system wide."
  exit 1
fi

SMASH_DIR_VAR=$1
EIGEN_DIR_VAR=$2
PYTHIA_DIR_VAR=$3

if [[ "${PYTHIA_DIR_VAR}" = '' ]]; then
    PYTHIA_OPTION_FOR_CMAKE=''
else
    PYTHIA_OPTION_FOR_CMAKE="-DPythia_CONFIG_EXECUTABLE=${PYTHIA_DIR_VAR}"
fi

export SMASH_DIR=$SMASH_DIR_VAR
EXAMPLE_DIR=${SMASH_DIR_VAR}/examples/using_SMASH_as_library

fail_and_rm_build() {
  echo "$*" >&2
  cd .. && rm -r build
  exit 1
}

cd $EXAMPLE_DIR
mkdir -p build && cd build
# Assume the default compiler is fine or that `CC` and `CXX` environment variables are set
cmake .. -DCMAKE_INSTALL_PREFIX=$EIGEN_DIR_VAR $PYTHIA_OPTION_FOR_CMAKE || fail_and_rm_build "Failed to setup cmake for SMASH libary examples"
make -j$(nproc) || fail_and_rm_build "Failed to build SMASH libary examples"

mkdir data  # create output directory for smash wrapper
./dynamic_fluid || fail_and_rm_build "Falhei"

#cd .. && rm -r build # Clean-up
