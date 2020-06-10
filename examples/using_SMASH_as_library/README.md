This is an example, how to set up a project using SMASH as a library.


## Prequisites

This example assumes that SMASH is installed and compiled, therefore, all the libraries it needs are already there.

## Creating the project

1. Set the SMASH_DIR environment variable by executing

        export SMASH_DIR=[...]/smash

- Copy the cmake files to find the libraries for the project

      export MY_PROJECT_DIR=[...]
      cp -r $SMASH_DIR/cmake $MY_PROJECT_DIR/

- Do the coding for your project (an example is given here) and run

      mkdir build && cd build
      cmake $MY_PROJECT_DIR -DCMAKE_INSTALL_PREFIX=[...]/eigen3 -DPythia_CONFIG_EXECUTABLE=[...]/pythia8302/bin/pythia8-config
      make
