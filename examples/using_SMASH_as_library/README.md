# Using SMASH as a library in a CMake project

This is meant as guidance to set up a project which uses SMASH as a library.
The examples included in this folder show how SMASH can be used as library to make use of specific functions of SMASH or how it can be wrapped as a whole.
The two main interface functions to be used when using SMASH as a library to setup and initialize are found and documented in _library.h_ in the SMASH source and are are used in the wrapper example.

## Prerequisites

This example assumes that SMASH is installed and, therefore, that all the libraries it needs are already there.
In particular, you can use `make install` from SMASH build folder to install SMASH, i.e. to copy where you wish all needed ingredients to use SMASH as a library.
Refer to the [INSTALL](../../INSTALL.md) file for more information about how to install SMASH and, possibly, customize the installation directory.

## Creating the CMake project that uses SMASH as a library

In order to set up your project, you need CMake to locate SMASH installation in order to then properly set include directories and linking libraries for your targets.
You can have a look to how this is done in the _CMakeLists.txt_ example file in this folder and get inspired.
However, it is crucial to be aware of the following needed steps.

1. Set the `SMASH_INSTALL_DIR` environment variable, which will be used by CMake to locate SMASH:
   ```bash
   export SMASH_INSTALL_DIR="${HOME}/.local"
   ```
   Optionally, you can use the `SMASH_VERSION` environment variable to pick the desired SMASH version, e.g. in case more than a version is installed.
   Its content should be the version number, e.g. `export SMASH_VERSION=3.0`.
   If you want to rely on a minimum or fixed version and make this a requirement of your project, add a version specification to your CMake code which finds SMASH.
   Use e.g.
   ```cmake
   find_package(SMASH 3.0 EXACT REQUIRED)
   ```
   to exclusively request version `SMASH-3.0` in your project (dropping `EXACT` would request _at least_ version `3.0`).
   Check out the documentation of [the `[version]` argument](https://cmake.org/cmake/help/latest/command/find_package.html#basic-signature) of the CMake `find_package` function for more information.

2. Your project will need some custom CMake files shipped by SMASH which serve to find external software on which SMASH relies.
   These are in the ***cmake*** folder and you can copy its content either to your project or to a different folder.
   For example, you could
   ```bash
   export MY_PROJECT_DIR='/path/to/your/project'
   cp -r /path/to/SMASH/cmake "${MY_PROJECT_DIR}"
   ```
   or, alternatively, e.g.
   ```bash
   mkdir -p ~/.cmake/modules
   cp /path/to/SMASH/cmake/* ~/.cmake/modules
   ```
   In both cases, you'll need to tell CMake where those files are.
   This can be done in your _CMakeLists.txt_ file via e.g.
   ```cmake
   list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")
   ```
   Alternatively, you can inform CMake about where the modules are at run time, e.g. via
   ```bash
   cmake -DCMAKE_MODULE_PATH=${HOME}/.cmak/modules ...
   ```
   or simply by permanently aliasing the `cmake` command, e.g. `alias cmake='cmake -DCMAKE_MODULE_PATH=~/.cmake/modules'`.

3. Setup your project using CMake in a standard way. For example:
   ```bash
   mkdir build && cd build
   cmake -DPythia_CONFIG_EXECUTABLE=path/to/pythia8316/bin/pythia8-config "${MY_PROJECT_DIR}"
   make
   ```

## Bash script in this folder

In this directory there is also a bash script that is internally used by SMASH for testing purposes.
It will install SMASH in a temporary folder and then build and run the examples present in this folder.
Therefore, if you read through the script, you will also find the instructions explained above.
However, be aware that you will also see some operations and checks on variables that are in general not needed; do not be confused by that and just refer to the operations above in an external project.
