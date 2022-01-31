## Introduction

This is an example of C++ code to read the SMASH HepMC3_treeroot output,
kindly provided by the HepMC3 developer Andrii Verbytskyi <andrii.verbytskyi@mpp.mpg.de>,
with minor additional modifications.

## Prerequisites

This example assumes that the user has successfully produced SMASH output files in HepMC3_treeroot format,
i.e. `SMASH_HepMC_particles.treeroot` and/or `SMASH_HepMC_collisions_conv.treeroot`,
and that ROOT is properly installed, so that the compiler can find the headers and the linker the shared libraries.
In the case of an installation from self-compiled sources, executing:
`source <path_to_ROOT_installation>/bin/thisroot.sh`
within the same shell instance that will be later used to try the example does usually the job.

## Compilation:

`g++ -o simple_treeroot_read_standalone.cc.{exe,cxx} $(root-config --cflags --glibs)`

## Execution:

Assuming that the output file is in the same directory as the newly created executable
(otherwise the appopriate path must be provided):

`./simple_treeroot_read_standalone.exe SMASH_HepMC_particles.treeroot`
