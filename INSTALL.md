# Install

---
* [More information about prerequisites](#prerequisites)
* [Customizing SMASH compilation](#building)
* [Customizing SMASH installation](#customizing)
* [Frequently asked questions](#faq)
   1. [SMASH does not compile. What should I do?](#not-compile)
   2. [SMASH crashes with 'illegal instruction'. Why?](#illegal-instruction)
   3. [SMASH output slightly differs on different platforms. Why?](#fuse-math-expr)
   4. [SMASH does not compile with pre-compiled ROOT binaries. What should I do?](#precompiled-root)
   5. [I run out of disk space compiling the code. Why?](#out-of-disk-space)
   6. [Can I disable tests or part of them?](#disabling-tests)
   7. [How can I use a different compiler?](#different-compilers)
   8. [How to use the LLVM implementation of the C++ standard library?](#llvm-STL)
   9. [How can I use SMASH as an external library?](#how-can-i-use-smash-as-an-external-library)
   10. [Can I disable ROOT or HepMC support?](#disable-root-hempc)
   11. [ROOT or HepMC are installed but CMake does not find them. What should I do?](#root-hepmc-not-found)

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

<a id="pythia"></a>

### Building Pythia

SMASH is tightly coupled to Pythia and thus requires a specific version, which is currently `8.316`.
If the required version is not already installed or if there are issues with the available one, it is recommended to build Pythia with similar flags as used for SMASH, like in the example below.

To download and build the needed version of Pythia, use the following commands:
```console
wget https://pythia.org/download/pythia83/pythia8316.tgz
tar xf pythia8316.tgz && rm pythia8316.tgz
cd pythia8316
./configure --cxx-common='-std=c++17 -march=native -O3 -fPIC -pthread'
make
```
To tell `cmake` where to find Pythia while building SMASH see the [**Building SMASH**](#building) section.

Note that although Pythia is statically linked into SMASH, access to `share/Pythia8/xmldoc` is required at runtime.

If you plan to build SMASH using the LLVM implementation of the standard C++ library, you should make sure that Pythia is built so as well, passing `-stdlib=libc++` together with the other flags to the `--cxx-common` option of the _configure_ script.
Of course, you should then use Clang compiler to build Pythia and this can either be achieved via the `CXX` environment variable (set via `export CXX=/path/to/clang++`) or by setting the variable for the Pythia configuration, e.g.
```console
CXX=clang++ ./configure --cxx-common='-std=c++17 -stdlib=libc++ -march=native -O3 -fPIC -pthread'
```

#### Remarks for Apple users

1. The `wget` command is not directly available on OSX.
   Although this can be easily installed e.g. via `brew install wget`, to download Pythia it is enough to use the `curl` command (see example below).

2. On recent Apple machines provided with M1 (ARM) processors, no `gcc` compiler is available and `clang` is to be used.
   The compiler flag `-march=native` is not defined and has to be dropped.

The commands to build Pythia on a M1 Apple machine become:
```console
curl https://pythia.org/download/pythia83/pythia8316.tgz -o pythia8316.tgz
tar xf pythia8316.tgz && rm pythia8316.tgz
cd pythia8316
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

This creates a folder named `gsl-[version_number]`.
In the following terminal commands, the absolute path to this folder is called `${GSL}`.
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

Clone the SMASH repository with the help of `git clone`.
The bare minimum commands needed to build SMASH from within the cloned repository are:
```console
mkdir build
cd build
cmake -DPythia_CONFIG_EXECUTABLE=[...]/pythia8316/bin/pythia8-config ..
make
```

However, CMake offers plenty of possible customizations, partly natively and partly created ad-hoc for SMASH.
In the following, the relevant explanation about these can be found and users should collect the relevant information for their case and build the appropriate `cmake` command according to their needs.

Please note that the `make` command builds everything (executables, tests, and libraries) and this might take a while.
You can use `make smash` if you are interested only in the executables or use `make smash_shared` to exclusively build the libraries (needed e.g. in another project [using SMASH as a library](#how-can-i-use-smash-as-an-external-library)).

#### Alternatives to specify the installation directory of Pythia

A few GNU/Linux distributions provide pre-built Pythia binaries without `pythia8-config` binary.
In this case, using the `-DPythia_CONFIG_EXECUTABLE` option as shown above is not possible and the top installation directory of Pythia containing `lib` has to be specified in either of the following ways:

*  Either set the bash environment variables `PYTHIA8` or `PYTHIA_ROOT_DIR` (e.g. `export PYTHIA_ROOT_DIR=/opt/pythia8316`) or
*  use the CMake `-DPYTHIA_ROOT_DIR` option (e.g. `cmake -DPYTHIA_ROOT_DIR=/opt/pythia8316 ..`).

If no variables are set and no options are passed, CMake searches for Pythia under the default path `/usr`.
To check which environment variables related to PYTHIA are currently set, use e.g. `printenv | grep PYTHIA`.


<a id="customizing"></a>

## SMASH installation

Installing SMASH gives the advantage that it is possible to simply use the `smash` command from anywhere in order to run SMASH (provided that the installation `bin` directory is contained in the `PATH` system environment variable).
Users planning to use SMASH as a library are encouraged to install SMASH, in order to disentangle the SMASH codebase source code from the version used as a library from other software.

The default installation done via `make install` installs SMASH into `/usr/local`.
If you want to change the installation directory, use the `CMAKE_INSTALL_PREFIX` CMake option to specify a new location.
For example, if you want to install SMASH in `~/.local`, use
```console
cmake -DCMAKE_INSTALL_PREFIX=${HOME}/.local ..
make install
```
You might need to use `sudo` for a system-wide installation.

By specifying `-DCMAKE_INSTALL_PREFIX=prefix`,
* `prefix/bin` will contain programs - e.g., `smash`,
* `prefix/lib/smash-X.Y-suffix` will contain libraries - e.g., `libsmash.so`,
* `prefix/include/smash-X.Y-suffix` will contain headers, and
* `prefix/share/smash-X.Y-suffix` will contain data files.

The string `X.Y-suffix` refers to the SMASH version and the `-suffix` to the Git commit at which the codebase was installed.
This won't appear for SMASH stable releases.
Such a suffix is also used for the executable, but a symbolic link to it named `smash` is also created in the binary folder, so that the `smash` command can be used after the installation.

On Unix OS you might be used to use [the `DESTDIR` mechanism](https://cmake.org/cmake/help/latest/envvar/DESTDIR.html) to relocate the whole installation, but refer to the CMake documentation to be sure that this is what you want.

**NOTE:** CMake does not always support multi-core builds using the installation target and, on some machines, `make -jN install` results in a warning similar to

> gmake[1]: warning: jobserver unavailable: using -j1.  Add '+' to parent make rule.

which, in turn, informs you that the build for the installation will be done with one single core.
It is likely that future CMake versions will fix this, but a simple work around is to first build SMASH and then install it: `make -jN smash && make install`.


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

### SMASH output slightly differs on different platforms. Why?

By default, SMASH is compiled optimizing for the native architecture, i.e. that where the compilation is done.
If hardware supports different optimizations on different machines, the default compilation of SMASH can indeed lead to (physically irrelevant) differences in the outcome of the same identical run.
The major (and possibly only) source of discrepancy that is worth mentioning is the usage of contracted mathematical expressions, which in SMASH default compilation are allowed, if the hardware supports them.
How to permit the compiler to form fused floating-point operations, such as fused multiply-add (FMA), depends on the compiler in use and you should check out which is the default value of the `-ffp-contract` flag in your compiler documentation.
Allowing contracted expression will produce slightly more precise results, since there is usually one rounding less operation.
However, exact (bit-wise) reproducibility of results is then not guaranteed and you will likely have (physically irrelevant) discrepancies across different platforms.

If this is not acceptable for your use case (as it is not for SMASH tests, where this feature is disabled), when setting up SMASH compilation, you can disable this feature, by passing the correct flag to CMake e.g. via
```console
cmake -DPythia_CONFIG_EXECUTABLE=[...] -DCMAKE_CXX_FLAGS="-ffp-contract=off" -DCMAKE_C_FLAGS="-ffp-contract=off"` ..
```

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

<a id="disabling-tests"></a>

### Can I disable tests or part of them?

SMASH is shipped with many tests of different type.
Most of them are unit and integration tests which have been kept separated from functional tests.
When CMake configure the project (i.e. when running `cmake` from the ***build*** folder), unit and integration tests are always setup to be later compiled and each of them has an executable associated.
However, [as mentioned](#out-of-disk-space), one is not obliged to compile everything, as it is possible to pass a given target to `make`.

On the other hand, functional tests are only setup if the CMake option `ENABLE_FUNCTIONAL_TESTS` is set to `ON`, which is **NOT** the default case.
This is due to the fact that, in order to setup these tests, CMake will create a Python virtual environment installing requirements in it and this might take some time (usually the first time only) if some required packages (e.g. `pandas`) need to be built from source.
If these tests are not enabled, they will not be visible and cannot be executed when running `ctest`.
Pass `-DENABLE_FUNCTIONAL_TESTS=ON` to `cmake` in order to include these tests when setting up the project.
Functional tests are written in Python and they run SMASH as black-box.
In this sense they do not need to be compiled, but the `smash` executable has to be created prior to their execution.

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

It is worth noticing that, when a GCC system-wide installation is present (e.g. under GNU/Linux Ubuntu), Clang compiler is by default locating it, in order to use its include files and its STL library implementation.
If several GCC installations are found, the most recent version is selected.
You can check which installation is selected by running `clang++ -v`.
It is important to be aware of this mechanism to possibly explain unexpected (but also unlikely) failures or simply different behavior.
For example, if the system default GCC compiler is not the latest among the available installations, Clang will use a different version of the STL from the GCC compiler.
Furthermore, if the most recent GCC compiler installation were lacking some needed file(s), Clang would simply fail either at compilation or at linking time.

**NOTE:** The FPE environment only works with GCC, so e.g. you won't get back-traces from floating point traps with Clang.

<a id="llvm-STL"></a>

### How to use the LLVM implementation of the C++ standard library?

In case the system default implementation of the C++ standard library is e.g. that shipped with GCC, this will still be used even when requesting CMake to use Clang as compiler.
However, it is possible to request to use the LLVM implementation using the CMake `CLANG_USE_LIBC++` option, e.g.
```console
CC=clang CXX=clang++ cmake -DPythia_CONFIG_EXECUTABLE=[...] -DCLANG_USE_LIBC++=ON ..
```
Still, when CMake performs some preliminary checks on the compiler, these are done using the default STL implementation found by Clang (i.e. that of the most recent GCC system-wide installation).
If you want to go 100% LLVM or in case some default file looked for by Clang is missing, you also need to set the corresponding CMake flag via
```console
CC=clang CXX=clang++ CXXFLAGS=-stdlib=libc++ cmake -DPythia_CONFIG_EXECUTABLE=[...] -DCLANG_USE_LIBC++=ON ..
```
which can also given as `-DCXXFLAGS="-stdlib=libc++"` argument to the `cmake` command.

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

**NOTE:** Remember to compile Pythia using LLVM implementation, too, as [previously described](#pythia).

<a id="how-can-i-use-smash-as-an-external-library"></a>

### How can I use SMASH as an external library?

The recommended way to use SMASH as a library in another software is to first install it according to the [instructions in the install section](#customizing).
This will build and copy all necessary SMASH files to the installation folder and, therefore, prevent unexpected surprises in your software behavior in case the SMASH source code changes (e.g. because of git operations).
In your software you can then use the files in the SMASH installation directory and manually pass them to the compiler.

Another possibility would be to build SMASH following the [instructions in the building section](#building).
Please be aware that in this case, using `make` (builds executables, tests, and libraries) or `make smash_shared` (exclusively builds the libraries) will both work, whereas `make smash` will not since this only builds the `smash` executable.
Be aware that this build is just local and you would have to refer to the path of it in the project that needs to access the SMASH libraries.

However, we encourage you to set up your project with CMake, too.
In that case, you can use the _FindSMASH.cmake_ module shipped in the ***cmake*** folder.
Refer to the *examples/using_SMASH_as_library/CMakeLists.txt* file for an example.
There are two important aspects to mention, in order to let the CMake `find_package(SMASH)` command succeed:
1. The folder where *FindSMASH.cmake* is needs to be in the `CMAKE_MODULE_PATH` CMake variable;
2. The `SMASH_INSTALL_DIR` environment variable must be correctly set.

If you are interested how to set up SMASH as an external library for your project, check out the [example how to do this in a CMake project](https://github.com/smash-transport/smash/tree/main/examples/using_SMASH_as_library).

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
