# Using SMASH as a library in a CMake project

This is meant as guidance to set up a project which uses SMASH as a library.
The examples included in this folder show how SMASH can be used as library to make use of specific functions of SMASH or how it can be wrapped as a whole.
The two main interface functions to be used when using SMASH as a library to setup and initialize are found and documented in _library.h_ in the SMASH source and are used in the wrapper example.

## Prerequisites

This example assumes that SMASH is installed and, therefore, that all the libraries it needs are already there.
In particular, you should use `make install` from SMASH build folder to install SMASH, i.e. to copy where you wish all needed ingredients to use SMASH as a library.
Refer to the [INSTALL](../../INSTALL.md) file for more information about how to install SMASH and, possibly, customize the installation directory.

With the modern CMake setup, SMASH installs CMake package configuration files and exposes components that can be requested by client projects.

## Locating the SMASH installation in your CMake project

To set up your project, you need CMake to locate the SMASH installation to then import and use the CMake targets offered by SMASH.
You can have a look to how this is done in the _CMakeLists.txt_ example file in this folder and get inspired.
However, it should be as simple as adding
```cmake
find_package(SMASH [<version>] [REQUIRED] [COMPONENTS <components>...])
```
to your CMakeLists file where the offered components are `Core`, `ROOT`, `HepMC3` and `Rivet`.
Each of those offers a target that gets imported if SMASH is correctly found.

### What if CMake complains about being unable to find SMASH?

First of all CMake will provide a descriptive hint.
Read it carefully and follow the advice given.
Most of the time, it happens because you might have installed SMASH in a non-standard location.
In such a case you need to help CMake locating SMASH, e.g. using the `CMAKE_PREFIX_PATH` environment variable or the `-DCMAKE_PREFIX_PATH` command line option to `cmake`.
For any `Package` CMake is also considering as a hint the `Package_DIR` folder, so you could specify the `SMASH_DIR` environment variable or pass `-DSMASH_DIR` to `cmake`.

### The imported targets

In the modern CMake spirit, after having called `find_package`, you should be ready to go.
SMASH components import the following targets:

| Component | Targets imported |
| :-------: | :--------------: |
| `Core`    | `SMASH::smash` |
| `ROOT`    | `SMASH::smash`, `SMASH::smash_root` |
| `HepMC3`  | `SMASH::smash`, `SMASH::smash_hepmc3` |
| `Rivet`   | `SMASH::smash`, `SMASH::smash_hepmc3`, `SMASH::smash_rivet` |

Any of these can be used to e.g. link your executable via
```cmake
target_link_libraries(my_executable PRIVATE SMASH::smash)
```
It is important to remark that all the non-core targets are an extensions of the core one.
Therefore, if you need ROOT functionalities in addition to the core SMASH library you can simply link as
```cmake
target_link_libraries(my_executable PRIVATE SMASH::smash_root)
```
and there is no need to specify there `SMASH::smash` as well.
Note that even though _SMASHConfig.cmake_ offers further targets, these are not guaranteed to be stable in time and might be changed (e.g. renamed) in the future.

### CMake variables imported when locating SMASH and SMASH dependencies

**After** a successful call to `find_package(SMASH ...)`, you will have some useful CMake variable available in your project.

* `SMASH_INPUT_FILES_DIR` provides the installation path to its input files.
  This can be used, for example, to pass the input directory to your code:
  ```cmake
  target_compile_definitions(project_requirements
                             INTERFACE SMASH_INPUT_DIR=\"${SMASH_INPUT_FILES_DIR}\")
  ```

* `Pythia_CONFIG_EXECUTABLE` provides the path to the Pythia configuration executable used to install SMASH.
  Since the SMASH installation is also installing CMake modules to locate dependencies that are not CMake-friendly (exactly like Pythia), it needs to keep track about how SMASH was able to locate them in first place.
  Therefore, some look-up variables are kept with the great advantage that **your project won't need to bother about this aspect**.
  The same applies to GSL, Eigen3 (mandatory dependencies) and ROOT, HepMC3 and Rivet (optional dependencies): They will be located exactly in the same way as they have been by SMASH.
  Finally, third-party libraries frozen inside SMASH (YAML, Cuba and Einhard) can be also be used in your project and linking to SMASH targets will correctly link them.

### Further functionality

You can pass to `find_package` any further command line option (refer to [CMake official documentation](https://cmake.org/cmake/help/latest/command/find_package.html)).
For example, if you want to rely on a minimum or fixed version and make this a requirement of your project, add a version specification to your CMake code which finds SMASH.
Use e.g.
```cmake
find_package(SMASH 3.4 EXACT REQUIRED)
```
to exclusively request version `SMASH-3.4` in your project (dropping `EXACT` would request _at least_ version `3.4`).
At the moment, specifying a version will make CMake accept any version larger than or equal to that specified.
Said differently, any SMASH version after the requested one will be considered compatible and hence accepted by CMake.


## Bash script in this folder

In this directory there is also a bash script that is internally used by SMASH for testing purposes.
It will install SMASH in a temporary folder and then build and run the examples present in this folder.
Therefore, if you read through the script, you will also find the instructions explained above.
However, be aware that you will also see some operations and checks on variables that are in general not needed; do not be confused by that and just refer to the operations above in an external project.
