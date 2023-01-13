# Install

---
* [More information about prerequisites](#prerequisites)
* [Customizing SMASH compilation](#building)
* [Customizing SMASH installation](#customizing)
* [Frequently asked questions](#faq)
   1. [SMASH does not compile. What should I do?](#not-compile)
   2. [SMASH crashes with 'illegal instruction'. Why?](#illegal-instruction)
   3. [I want the compiler to use contracted floating-point operations. What should I do?](#fuse-math-expr)
   4. [SMASH does not compile with pre-compiled ROOT binaries. What should I do?](#precompiled-root)
   5. [I run out of disk space compiling the code. Why?](#out-of-disk-space)
   6. [How can I use a different compiler?](#different-compilers)
   7. [How to use the LLVM implementation of the C++ standard library?](#llvm-STL)
   8. [Can I disable ROOT or HepMC support?](#disable-root-hempc)
   9. [ROOT or HepMC are installed but CMake does not find them. What should I do?](#root-hepmc-not-found)

---

## Preliminary note

In the following, whenever the `make` command is used, you can specify the `-j` option followed by the number of CPU cores that should be used in parallel, e.g. `-j3` to use three cores.
This will speed up compilation considerably.



<a id="prerequisites"></a>

## About prerequisites

All UNIX-like operating systems offer package managers and it is highly encouraged to use these whenever possible.
To mention a couple, you might use [Homebrew](https://brew.sh) on Apple machines or [APT](https://wiki.debian.org/Apt) on Debian/Ubuntu distributions.
The installation of required compilers, CMake, GSL and Eigen3 libraries (just to mention some) should be straightforward with a package manager.

It is very likely that a compiler among those supported is available out of the box on your operating system.
However, some of the SMASH prerequisites are less likely to be already available, Pythia first of all.

### Building Pythia

SMASH is tightly coupled to Pythia and thus requires a specific version, which is currently `8.307`.
If the required version is not already installed or if there are issues with the available one, it is recommended to build Pythia with similar flags as used for SMASH, like in the example below.

To download and build the needed version of Pythia, use the following commands:
```console
wget https://pythia.org/download/pythia83/pythia8307.tgz
tar xf pythia8307.tgz && rm pythia8307.tgz
cd pythia8307
./configure --cxx-common='-std=c++17 -march=native -O3 -fPIC -pthread'
make
```
To tell `cmake` where to find Pythia while building SMASH see the [**Building SMASH**](#building) section.

Note that although Pythia is statically linked into SMASH, access to `share/Pythia8/xmldoc` is required at runtime.

If you plan to build SMASH using the LLVM implementation of the standard C++ library, you should make sure that Pythia is built so as well, passing `-stdlib=libc++` together with the other flags to the `--cxx-common` option of the _configure_ script.

#### Remarks for Apple users

1. The `wget` command is not directly available on OSX.
   Although this can be easily installed e.g. via `brew install wget`, to download Pythia it is enough to use the `curl` command (see example below).

2. On recent Apple machines provided with M1 (ARM) processors, no `gcc` compiler is available and `clang` is to be used.
   The compiler flag `-march=native` is not defined and has to be dropped.

The commands to build Pythia on a M1 Apple machine become:
```console
curl https://pythia.org/download/pythia83/pythia8307.tgz -o pythia8307.tgz
tar xf pythia8307.tgz && rm pythia8307.tgz
cd pythia8307
./configure --cxx-common='-std=c++17 -O3 -fPIC -pthread'
make
```


### Installing Eigen

Usually it is possible to install Eigen with a package manager (it requires admin privileges) and in this case CMake should be able to find the header files without the need of any additional option.
For example, on an Apple machine you have the possibility to install Eigen via `brew install eigen`, while, under GNU/Linux Ubuntu, via `sudo apt-get install libeigen3-dev`.

If for some reason this is not a viable approach, then you can still proceed to a manual installation as described in the following.

#### Getting Eigen header files into a custom location

Let's assume Eigen headers will be unpacked in `${HOME}`.
Download the latest stable release of `Eigen` from [the official website](http://eigen.tuxfamily.org) and unpack it via
```console
tar -xf "[latest-eigen].tar.gz" -C "${HOME}"
```
To tell CMake where to find Eigen header files while building SMASH, pass the path to them adding the option `-DCMAKE_PREFIX_PATH=$HOME/[latest-eigen]/` to the `cmake` command in the SMASH setup (cf. [Building SMASH](#building) section).


### Installing and enabling Rivet support

[Rivet website](https://rivet.hepforge.org/) is pretty complete and general information can be found there.
Since the interface with SMASH has been tested with version 3.1.4, if possible, this version should be installed and used.

An [installation script](https://gitlab.com/hepcedar/rivet/-/blob/release-3-1-x/doc/tutorials/installation.md) exists and can be downloaded e.g. via
```console
wget https://gitlab.com/hepcedar/rivetbootstrap/raw/3.1.4/rivet-bootstrap
```
It provides a convenient way to install Rivet and its dependencies (HepMC3 is among those, but, if you have already installed it, you can edit the script so that Rivet uses your installation).

Please note that, every time Rivet is used, some environment variables must be set in advance.
The `rivetenv.sh` script, in the Rivet installation directory, takes care of this step, if sourced
```console
source "[...]/rivetenv.sh"
```
where `[...]` is not a command, but a shorthand for the path of the directory in which Rivet is installed.

If Rivet (with all its dependencies) is installed and the environment variables are set, it will be automatically detected by CMake during SMASH setup.


### Making and using a custom GSL installation

Download and unpack the latest GSL version via e.g.
```console
wget ftp://ftp.gnu.org/gnu/gsl/gsl-latest.tar.gz
tar -zxvf gsl-latest.tar.gz
```

This creates a folder named `gsl-[version_number]`, which is called `${GSL}` here.
To build and install the downloaded source use
```console
cd "${GSL}"
./configure --prefix "${GSL}"
make
make install
```

When setting up SMASH, use the `-DGSL_ROOT_DIR=${GSL}` CMake option to let CMake find your custom GSL installation.



<a id="building"></a>

## Building SMASH

The bare minimum needed to build SMASH from within its repository reads:
```console
mkdir build
cd build
cmake -DPythia_CONFIG_EXECUTABLE=[...]/pythia8307/bin/pythia8-config ..
make
```

However, CMake offers plenty of possible customizations, partly natively and partly created ad-hoc for SMASH.
In the following, the relevant explanation about these can be found and users should collect the relevant information for their case and build the appropriate `cmake` command according to their needs.

#### Alternatives to specify the installation directory of Pythia

A few GNU/Linux distributions provide pre-built Pythia binaries without `pythia8-config` binary.
In this case, using the `-DPythia_CONFIG_EXECUTABLE` option as shown above is not possible and the top installation directory of Pythia containing `lib` has to be specified in either of the following ways:

*  Either set the bash environment variables `PYTHIA8` or `PYTHIA_ROOT_DIR` (e.g. `export PYTHIA_ROOT_DIR=/opt/pythia8307`) or
*  use the CMake `-DPYTHIA_ROOT_DIR` option (e.g. `cmake -DPYTHIA_ROOT_DIR=/opt/pythia8307 ..`).

If no variables are set and no options are passed, CMake searches for Pythia under the default path `/usr`.
To check which environment variables related to PYTHIA are currently set, use e.g. `printenv | grep PYTHIA`.


<a id="customizing"></a>

## SMASH installation

Installing SMASH gives the advantage that it is possible to simply use the `smash` command from anywhere in order to run SMASH.

The default installation done via `make install` installs SMASH into `/usr/local`.
If you want to change the installation directory, use the `CMAKE_INSTALL_PREFIX` CMake option to specify a new location.
For example, if you want to install SMASH in `~/.local`, use
```console
cmake -DCMAKE_INSTALL_PREFIX=${HOME}/.local ..
make install
```

By specifying `-DCMAKE_INSTALL_PREFIX=prefix`,
* `prefix/bin` will contain programs - e.g., `smash`,
* `prefix/lib` will contain libraries - e.g., `libsmash.so`,
* `prefix/include/smash` will contain headers, and
* `prefix/share/smash` will contain data files.



<a id="faq"></a>

## FAQ

If you do not find help in the following answers, feel free to browse [closed issues](https://github.com/smash-transport/smash/issues?q=is%3Aissue+is%3Aclosed) or open a new one [on GitHub](https://github.com/smash-transport/smash/issues).

<a id="not-compile"></a>

### SMASH does not compile. What should I do?

If compilation fails (especially after changing a library), using a fresh build folder can sometimes fix the problem.

<a id="illegal-instruction"></a>

### SMASH crashes with "illegal instruction". Why?

If running SMASH fails with "illegal instruction", it is likely due to running SMASH on a different platform than the one where it was compiled.
By default, SMASH is compiled with platform-specific optimizations, implying that the binary only works on such platforms.

There are various possible ways to fix this issue:
1. Make sure to compile SMASH on the platform where you run it.
   When running SMASH on a computing cluster, this may require making the compilation step part of the job script, if e.g. the node architecture is different from that of the login or submission node.
2. Only run SMASH on platforms similar to the one where it was compiled.
   When running SMASH on a computing cluster, this may require restricting the jobs to the correct platform.
3. Find out for which platform you need to compile SMASH and specify it setting up the compilation.
   For this purpose, use CMake command line options, as e.g. `-DCMAKE_CXX_FLAGS="-march=x86-64" -DCMAKE_C_FLAGS="-march=x86-64"`.
4. As last resort, compile SMASH without machine-specific optimizations by removing `-march=native` from `CMakeLists.txt`.
   This is the easiest solution, however it results in a less efficient executable.
   Note that the same applies to any other libraries you compile with `-march=native`, for instance Pythia.

<a id="fuse-math-expr"></a>

### I want the compiler to use contracted floating-point operations. What should I do?

How to permit the compiler to form fused floating-point operations, such as fused multiply-add (FMA), depends on the compiler you use and you should check out which is the correct value of the flag `-ffp-contract` flag in the compiler documentation.
Allowing contracted expression will produce slightly more precise results, since there is usually one rounding less operation.
However, exact (bit-wise) reproducibility of results is then not guaranteed and you will likely have (physically irrelevant) discrepancies across different platforms.

If this is acceptable for your use case, when setting up SMASH compilation, pass the correct flag to CMake.
For LLVM `clang` compiler, use e.g. `-DCMAKE_CXX_FLAGS="-ffp-contract=on" -DCMAKE_C_FLAGS="-ffp-contract=on"`.
The GNU compiler (checked up to `g++-12.2`) implements the `-ffp-contract=on` flag as `-ffp-contract=off` and the way to allow fused operations is to restore the `-ffp-contract=fast` default value.

<a id="precompiled-root"></a>

### SMASH does not compile with pre-compiled ROOT binaries. What should I do?

If compilation of SMASH in combination with a pre-compiled ROOT binary fails, please install and compile ROOT locally from source (see http://root.cern.ch) and compile SMASH again in a clean build directory.

<a id="out-of-disk-space"></a>

### I run out of disk space compiling the code. Why?

Please note that after compilation the `smash` directory (including `build`) might have a size of some GB.
By default, the unit tests are always compiled, which requires a lot of the disk space.
If disk space is restricted, consider to just run `make smash`, which will only compile the SMASH binary.
However, it is still recommended to run the unit tests at least once when compiling in a new environment to ensure that everything works as expected.
To see how to run the tests, see [CONTRIBUTING](CONTRIBUTING.md).

<a id="different-compilers"></a>

### How can I use a different compiler?

In order to use a particular compiler, e.g. Clang, you can permanently set the following environment variables
```console
export CC=clang
export CXX=clang++
```
or simply set them for the `cmake` command only via
```console
CC=clang CXX=clang++ cmake ..
```
Alternatively, the compiler can also be specified using CMake dedicated options,
```console
cmake -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ ..
```
**NOTE:** The FPE environment only works with GCC, so e.g. you won't get back-traces from floating point traps with Clang.

<a id="llvm-STL"></a>

### How to use the LLVM implementation of the C++ standard library?

In case the system default implementation of the C++ standard library is e.g. that shipped with GCC, this will still be used even when requesting CMake to use Clang as compiler.
However, it is possible to request to use the LLVM implementation using the CMake `CLANG_USE_LIBC++` option, e.g.
```console
CC=clang CXX=clang++ cmake -DPythia_CONFIG_EXECUTABLE=[...] -DCLANG_USE_LIBC++=ON ..
```
If the installation of the LLVM implementation is not in a standard place, you either need to set and export your `LD_LIBRARY_PATH` environment variable to the correct value, e.g.
```console
export LD_LIBRARY_PATH=/path/to/clang/installation/lib
```
or pass to the `cmake` command the option
```
-DCMAKE_EXE_LINKER_FLAGS="-Wl,-rpath -Wl,/path/to/clang/installation/lib"
```
where, of course, the path to clang installation must be a valid path.
All of this is needed to let the executable find the library ABI at run time.

<a id="disable-root-hempc"></a>

### Can I disable ROOT or HepMC support?

By default, if SMASH setup finds ROOT and/or HepMC, it will use them.
However, this feature can be disabled by using the CMake options `-DTRY_USE_ROOT=OFF` and/or `-DTRY_USE_HEPMC=OFF`.
For example,
```console
cmake -DTRY_USE_ROOT=OFF -DTRY_USE_HEPMC=OFF <source_dir>
```
will setup SMASH without ROOT and without HepMC support.

<a id="root-hepmc-not-found"></a>

### ROOT or HepMC are installed but CMake does not find them. What should I do?

If the ROOT or HepMC installation is not found, provide a hint to CMake about the install destination with
```console
cmake -DCMAKE_PREFIX_PATH=/path/to/root/or/HepMC/installation ..
```
Note that if multiple `CMAKE_PREFIX_PATH`s are necessary, a semicolon-separated list of directories can be specified.
