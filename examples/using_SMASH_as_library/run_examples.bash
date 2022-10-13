#!/bin/bash

# Eventually input
SMASH_DIR_VAR="~/smash-devel"
EIGEN_DIR_VAT="~/code/bin/eigen-eigen-b3f3d4950030/"
PYTHIA_DIR_VAR="~/code/bin/pythia8307/bin/pythia8-config"

export SMASH_DIR=$SMASH_DIR_VAR

mkdir -p build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=$EIGEN_DIR_VAR -DPythia_CONFIG_EXECUTABLE=$PYTHIA_DIR_VAR
make

./example
./example_rate_equations
mkdir -p data  # create output directory for smash wrapper
./example_smash_wrapper
