# SMASH README

This is the repository for the development of the SMASH (Simulating Many
Accelerated Strongly-interacting Hadrons) transport approach for the dynamical
description of heavy ion reactions. See [Phys. Rev. C 94, 054905
(2016)](https://arxiv.org/abs/1606.06642) for details.

See [CONTRIBUTING](CONTRIBUTING.md) for development hints. A complete User
Guide is found [here](https://fias.uni-frankfurt.de/~smash/extra/user/).

If Pythia is used, please cite the following references:

* T. Sjöstrand, S. Mrenna and P. Skands, JHEP05 (2006) 026,
  Comput. Phys. Comm. 178 (2008) 852.

## How to build SMASH

### Prerequisites

SMASH is known to compile and work with one of these compilers (which have the
required C++11 features):
- gcc >= 4.8
- clang >= 3.2

It requires the following tools & libraries:
- cmake >= 2.8.11
- the GNU Scientific Library >= 1.15
- the Eigen3 library for linear algebra (see http://eigen.tuxfamily.org)
- boost filesystem >= 1.49

Support for ROOT output is automatically enabled if a suitable version of ROOT
(>= 5.34) is found on the system.


### Building SMASH

Use the following commands too build SMASH in a separate directory:

    mkdir build
    cd build
    cmake ..
    make

To build in parallel on N cores:

    make -jN

To run it with specific settings:

    vi config.yaml
    ./smash


### Size of the code

Please note that after compilation the `smash` directory (including `build`)
has a size of about 4GB. If disk space is restricted, consider to just run

    make smash

which will only compile the SMASH binary. By default, the unit tests are also
compiled which requires a lot of the disk space. It is still recommended to run
the unit tests at least once when compiling in a new environment to ensure that
everything works as expected. To see how to run the tests see
[CONTRIBUTING](CONTRIBUTING.md).


### Changing the compiler

In order to use a particular compiler, you can set the following environment
variables:

    export CC=gcc
    export CXX=g++

Alternatively the compiler can also be specified to cmake like this:

    cmake .. -DCMAKE_CXX_COMPILER=g++


### Using the official clang version on macOS

The clang compiler version shipped with XCode or the Command Line Tools from
Apple in the past (before version 8.1) did not always support all C++ features
that were used. If problems with compiling on macOS arise, try the official
clang version. This can be done for example with [Homebrew](http://brew.sh):

    brew install llvm

After that you have to instruct cmake that it should use the newly installed
clang compiler and also tell the linker where the libraries are by altering some
cmake flags. Enter the cmake command like below:

    cmake .. -DCMAKE_CXX_COMPILER=/usr/local/opt/llvm/bin/clang++ -DCMAKE_EXE_LINKER_FLAGS="-L/usr/local/opt/llvm/lib -lc++abi" -DCMAKE_CXX_FLAGS=-I/usr/local/opt/llvm/include

Note: The FPE environment only works with gcc, so e.g. you won't get backtraces
from floating point traps with clang in general.


### Disabling ROOT support

Producing ROOT output requires ROOT installed (see http://root.cern.ch).
If ROOT is found, the support for ROOT output is automatically enabled.
In order to disable it, one can do the follwoing:

    cmake -DUSE_ROOT=OFF <source_dir>
    make


### Including Eigen header files from custom location

Let's assume Eigen headers will be unpacked in `$HOME`.

1. Download latest package `[latest-eigen].tar.gz` from http://eigen.tuxfamily.org

2. unpack: `tar -xf [latest-eigen].tar.gz -C $HOME`

3. in `smash/build/`, create build files with `cmake -DCMAKE_INSTALL_PREFIX=$HOME/[latest-eigen]/ ..`


### Using a custom GSL build

Download and unpack GSL:

    wget ftp://ftp.gnu.org/gnu/gsl/gsl-latest.tar.gz
    tar -zxvf gsl-latest.tar.gz

This creates a folder named `gsl-[version_number]` called `$GSL` here.

    cd $GSL
    ./configure --prefix $GSL
    make -jN
    make install

Here N is the number of cores to use in the `make` command. When compiling
SMASH, run `cmake` with

    cmake -DGSL_ROOT_DIR=$GSL ..

Note: In case of problems, make sure to start with a clean build folder.


## Running SMASH with example input files

SMASH ships with example configuration files for the collider, box, sphere and
list modus. By default, i.e. by running `./smash`, the simulation is set up by
means of a collider configuration file, called `config.yaml` and the default
particles and decaymodes files, `particles.txt` and `decaymodes.txt`. They are
located in `/input`.

Additionally, example configuration files for the box, sphere and list modus can
be found in the respective directories `/input/{box,sphere,list}`. In case of a
box simulation, the default particles and decaymodes files need to be modified
to allow for equilibration. These are also stored in `/input/box`. For the list
modus, an input list file to be read in is required. This file, `example_list0`,
is located in `/input/list`.

To run SMASH with a manually specified configuration file, use the `-i` command.
For example, for the sphere or list example file:

    ./smash -i ../input/sphere/config.yaml
    ./smash -i ../input/list/config.yaml

To further use non-default particles and decaymodes files, the `-p`
and `-d` options are necessary. For the default box, this means:

    ./smash -i ../input/box/config.yaml -p ../input/box/particles.txt -d ../input/box/decaymodes.txt

All command line options can be viewed with

    ./smash -h


## License

SMASH is licensed under the terms of the GNU General Public License, Version 3
or above. The build scripts are licensed under terms of the BSD 3-clause
license. See [LICENSE](LICENSE).

### Projects using SMASH

SMASH source and documentation are provided to check and
reproduce published results of the authors. Cooperation and joint projects with outside
researchers are encouraged and comparison to results by experimental collaborations
is supported. If you are interested in starting a project, please contact us to avoid
interference with current thesis topics. If your project involves changes to the code,
please refer to [CONTRIBUTING](CONTRIBUTING.md) for coding guidelines and
helpful tools.
