This is an example, how to set up a project using SMASH as a library.


## Prequisites

This example assumes that SMASH is installed and compiled, therefore, all the libraries it needs are already there.

## Creating project

- Set SMASH_DIR environment variable

  export SMASH_DIR=[...]/smash

- Copy cmake files for finding libraries to the project

  export MY_PROJECT_DIR=[...]
  cp -r $SMASH_DIR/cmake $MY_PROJECT_DIR/

- Do the coding for your project (an example is given here) and run

  mkdir build && cd build
  cmake $MY_PROJECT_DIR -DPythia_CONFIG_EXECUTABLE=[...]/pythia8230/bin/pythia8-config
  make
