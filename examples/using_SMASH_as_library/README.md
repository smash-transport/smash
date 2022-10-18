This is an example, how to set up a project using SMASH as a library. The two main interface functions
to be used when using SMASH as a library to setup and initialize are found and documented in library.h
in the SMASH source.

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
      cmake $MY_PROJECT_DIR -DCMAKE_INSTALL_PREFIX=[...]/eigen3 -DPythia_CONFIG_EXECUTABLE=[...]/pythia8307/bin/pythia8-config
      make


## Script

There is also a script found in this directory (run_examples_for_testing.bash) that is used for building and running the examples for testing. If you read through the script, you will also find the instructions explained above. In the script the build directory is automatically removed after running, but you can uncomment the corresponding lines, if you want to run it for yourself.
