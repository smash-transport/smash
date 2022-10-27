#!/bin/bash

if [[ $# -ne 3 ]]; then
  echo "Script requires exactly three arguments. Usage: $0 SMASH_DIR EIGEN_DIR PYTHIA_DIR_VAR"
  exit 1
fi

SMASH_DIR_VAR=$1
EIGEN_DIR_VAR=$2
PYTHIA_DIR_VAR=$3

export SMASH_DIR=$SMASH_DIR_VAR
EXAMPLE_DIR=${SMASH_DIR_VAR}/examples/using_SMASH_as_library

fail_and_rm_build() {
  echo "$*" >&2
  cd .. && rm -r build
  exit 1
}

cd $EXAMPLE_DIR
mkdir -p build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=$EIGEN_DIR_VAR -DPythia_CONFIG_EXECUTABLE=$PYTHIA_DIR_VAR -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ || fail_and_rm_build "Failed to setup cmake for SMASH libary examples"
make || fail_and_rm_build "Failed to build SMASH libary examples"

./example || fail_and_rm_build "Failed to run SMASH libary examples"
./example_rate_equations || fail_and_rm_build "Failed to run SMASH rate equation libary examples"
mkdir data  # create output directory for smash wrapper
./example_smash_wrapper || fail_and_rm_build "Failed to execute SMASH wrapper libary example"

cd .. && rm -r build # Clean-up
