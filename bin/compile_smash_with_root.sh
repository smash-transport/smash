#!/bin/bash

# Bash script to compile and run SMASH in combination with ROOT and HepMC3 on
# the gsi kronos cluster. ROOT is available in
# /cvmfs/it.gsi.de/modulefiles/ and needs only be loaded.
# To use the HepMC3 module, follow the instructions in the SMASH wiki. 

# For compilation, the paths to the SMASH source directory, to the eigen library
# as well as to the gsl and pythia directories are necessary and need to be
# passed from the command line.

if [[ "$#" -ne 4 ]]; then
    echo "Submit SLURM job for given target."
    echo
    echo "Usage: $0 SMASH_SRC_DIR EIGEN_SRC_DIR GSL_DIR PYTHIA_DIR"
    exit 1
fi

# exit when any command fails
set -e

# load ROOT
module use /cvmfs/it.gsi.de/modulefiles/
module load root/v6.06-06

# load HepMC3
module load use.own
module load hepmc3

smash_dir=$1
eigen_dir=$2
gsl_dir=$3
pythia_dir=$4

build_dir=$smash_dir/build_with_ROOT_and_HepMC

# Compile SMASH and run tests
cd $smash_dir \
&& cmake -DCMAKE_BUILD_TYPE=Release -DGSL_ROOT_DIR=$gsl_dir -DPythia_CONFIG_EXECUTABLE=$pythia_dir/bin/pythia8-config -DCMAKE_INSTALL_PREFIX=$eigen_dir -B$build_dir -H$smash_dir \
&& cd $build_dir \
&& make -j$SLURM_CPUS_ON_NODE

# Run tests
ctest -j$SLURM_CPUS_ON_NODE

# Run SMASH setup with ROOT and HepMC3 output enabled
output_dir=$build_dir/output_ROOT_HepMC
./smash -c "Output: {Particles: {Format: ["Root"]} }" -c "Output: {HepMC: {Format: ["ASCII"]} }" -o $output_dir

file_root=$output_dir/Particles.root
file_hepmc=$output_dir/SMASH_HepMC.asciiv3

# Verify that output files were created and SMASH did not crash
if [[ -f "$file_root" ]] && [[ -f "$file_hepmc" ]] && ! [[ -f "$output_dir/smash.lock" ]]; then
    printf '\n \u2714 SMASH ran successfully with Root and HepMC output enabled.\n\n'
else
    exit 1
fi
