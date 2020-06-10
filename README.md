# SMASH

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3484711.svg)](https://doi.org/10.5281/zenodo.3484711)

SMASH (Simulating Many Accelerated Strongly-interacting Hadrons) is a
relativistic hadronic transport approach for the dynamical description of heavy
ion reactions. Please see [Phys. Rev. C 94, 054905
(2016)](https://arxiv.org/abs/1606.06642) for details and cite this reference
together with the [software DOI](https://doi.org/10.5281/zenodo.3484711) for
the specific code version employed, if you are using SMASH. A BibTeX entry for
the software DOI is found on the respective Zenodo pages.


See [CONTRIBUTING](CONTRIBUTING.md) for development hints. A complete User Guide
can be found [here](https://theory.gsi.de/~smash/userguide/current/). For a more
detailed development documentation, see
[here](http://theory.gsi.de/~smash/doc/current/).


If Pythia is used, please cite the following references:
T. Sjöstrand, S. Mrenna and P. Skands,
[JHEP05 (2006) 026](https://arxiv.org/abs/hep-ph/0603175),
[Comput. Phys. Comm. 178 (2008)](https://arxiv.org/abs/0710.3820).


Report issues at https://github.com/smash-transport/smash/issues or contact us by email
at elfner@itp.uni-frankfurt.de.

## How to Build SMASH

### Prerequisites

SMASH is known to compile and work with one of these compilers (which have the
required C++11 features):
- gcc >= 4.8
- clang >= 3.2

It requires the following tools & libraries:
- cmake >= 2.8.11
- the GNU Scientific Library >= 2.0
- the Eigen3 library for linear algebra (see http://eigen.tuxfamily.org)
- boost filesystem >= 1.49
- Pythia = 8.302

Support for ROOT output is automatically enabled if a suitable version of ROOT
(>= 5.34) is found on the system.

### Building Pythia

SMASH is tightly coupled to Pythia and thus requires a specific version. Using
a different version than specified above may or may not work. It is recommended
to build Pythia with similar flags as used for SMASH:

    wget http://home.thep.lu.se/~torbjorn/pythia8/pythia8302.tgz
    tar xf pythia8302.tgz && rm pythia8302.tgz
    cd pythia8302
    ./configure --cxx-common='-std=c++11 -march=native -mfpmath=sse -O3 -fPIC'
    make

To tell cmake where to find Pythia, pass the path to the pythia8-config
executable as shown in the next section.

Note that although Pythia is statically linked into SMASH, access to
`share/Pythia8/xmldoc` is required at runtime.

### Including Eigen Header Files from Custom Location

Let's assume Eigen headers will be unpacked in `$HOME`.

1. Download latest package from http://eigen.tuxfamily.org

       [latest-eigen].tar.gz
       tar -xf [latest-eigen].tar.gz -C $HOME`

2. in `smash/build/`, create build files with `cmake -DCMAKE_INSTALL_PREFIX=$HOME/[latest-eigen]/ ..`


### Building SMASH

Use the following commands to build SMASH in a separate directory:

    mkdir build
    cd build
    cmake .. -DPythia_CONFIG_EXECUTABLE=[...]/pythia8302/bin/pythia8-config
    make

To build in parallel on N cores:

    make -jN

To run it with specific settings:

    vi config.yaml
    ./smash


### Troubleshooting

#### SMASH does not compile

If compilation fails (especially after changing a library), using a fresh build
folder can sometimes fix the problem.


#### SMASH crashes with "illegal instruction"

If running SMASH fails with "illegal instruction", it is likely due to running
SMASH on a different platform than the one where it was compiled. By default,
SMASH is compiled with platform-specific optimizations, implying that the binary
only works on such platforms.

There are three possible ways to fix this issue:

1. Make sure to compile SMASH on the platform where you run it. When running
   SMASH on a computing cluster, this requires making the compilation step part
   of the job script.
2. Only run SMASH on platforms similar to the one where it was compiled. When
   running SMASH on a computing cluster, this requires restricting the jobs to
   the correct platform.
3. Compile SMASH without machine-specific optimizations by removing
   `-march=native` from `CMakeLists.txt`. This is the easiest solution, however
   it results in less efficient code.

Note that the same applies to any other libraries you compile with
`-march=native`, like for instance Pythia.


### Size of the Code

Please note that after compilation the `smash` directory (including `build`)
has a size of about 4GB. By default, the unit tests are always compiled, which
requires a lot of the disk space. If disk space is restricted, consider to just run

    make smash

which will only compile the SMASH binary. It is still recommended to run
the unit tests at least once, when compiling in a new environment to ensure that
everything works as expected. To see how to run the tests, see
[CONTRIBUTING](CONTRIBUTING.md).


### Changing the Compiler

In order to use a particular compiler, you can set the following environment
variables:

    export CC=gcc
    export CXX=g++

Alternatively the compiler can also be specified to cmake like this:

    cmake .. -DCMAKE_CXX_COMPILER=g++

Note: The FPE environment only works with gcc, so e.g. you won't get backtraces
from floating point traps with clang.


### Disabling ROOT Support

Producing ROOT output requires ROOT installed (see http://root.cern.ch).
If ROOT is found, the support for ROOT output is automatically enabled.
In order to disable it, one can do the following:

    cmake -DUSE_ROOT=OFF <source_dir>
    make


### Using a Custom GSL Build

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

## Running SMASH with Example Input Files

SMASH ships with example configuration files for the collider, box, sphere, and
list modus. By default, i.e. by running `./smash`, the simulation is set up from
the collider configuration file, called `config.yaml`, and the default
particles and decaymodes files, `particles.txt` and `decaymodes.txt`. They are
located in `/input`.

Additionally, example configuration files for the box, sphere and list modus can
be found in the respective directories `/input/{box,sphere,list}`. In case of the
box simulation, in order to allow for equilibration, different default particles
and decaymodes files need to be used. These files are also provided in
`/input/box`. For the list modus, an input list file to be read in is required.
This file, `example_list0`, is located in `/input/list`.

To run SMASH with a non-default configuration file, use the `-i` command.
For example, for the sphere or list example file:

    ./smash -i ../input/sphere/config.yaml
    ./smash -i ../input/list/config.yaml

Further, using non-default particles and decaymodes files is necessary, and these
can be specified through the `-p` and `-d` options. In the box and the dileptons
example,
this means:

    ./smash -i ../input/box/config.yaml -p ../input/box/particles.txt -d ../input/box/decaymodes.txt
    ./smash -i ../input/dileptons/config.yaml -d ../input/dileptons/decaymodes.txt

All command line options can be viewed with

    ./smash -h

To run SMASH completely silently for production runs, we recommend to redirect
stdout to /dev/null, warnings and error messages will still be displayed.

    ./smash > /dev/null


## License

SMASH is licensed under the terms of the GNU General Public License, Version 3
or above. The build scripts are licensed under terms of the BSD 3-clause
license. See [LICENSE](LICENSE).

### Projects Using SMASH

SMASH source and documentation are provided to check and
reproduce published results of the authors. Cooperation and joint projects with
outside researchers are encouraged and comparison to results by experimental
collaborations is supported. SMASH can be used as a 3rd party library, for
examples see ${SMASH_DIR}/examples/ folder. If you are interested in starting a
project, please contact us to avoid interference with current thesis topics.
If your project involves changes to the code, please refer to
[CONTRIBUTING](CONTRIBUTING.md) for coding guidelines and helpful tools.
