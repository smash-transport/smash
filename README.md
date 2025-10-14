# SMASH

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3484711.svg)](https://doi.org/10.5281/zenodo.3484711)

SMASH (Simulating Many Accelerated Strongly-interacting Hadrons) is a relativistic hadronic transport approach for the dynamical description of heavy-ion reactions.
Please see [Phys. Rev. C 94, 054905 (2016)](https://arxiv.org/abs/1606.06642) for details and, if you are using SMASH, cite this reference together with the [software DOI](https://doi.org/10.5281/zenodo.3484711) for the specific code version employed.
A BibTeX entry for the software DOI is found on the respective Zenodo pages.

See [CONTRIBUTING](CONTRIBUTING.md) for development hints.
A complete [User Guide](https://theory.gsi.de/~smash/userguide/current/) as well as a more detailed [development documentation](http://theory.gsi.de/~smash/doc/current/) are available for the latest version of the code.
For documentation of older versions, refer to links in the [releases pages](https://github.com/smash-transport/smash/releases).

If Pythia is used, please cite the following references (both article and the codebase release you used):
* [_A comprehensive guide to the physics and usage of PYTHIA 8.3_](https://scipost.org/SciPostPhysCodeb.8), C. Bierlich et al; SciPost Phys. Codebases 8 (2022), DOI: `10.21468/SciPostPhysCodeb.8`, also available on [arXiv](https://arxiv.org/abs/2203.11601);
* [SciPost Phys. Codebases 8-r8.3](https://scipost.org/SciPostPhysCodeb.8-r8.3) (2022), DOI: `10.21468/SciPostPhysCodeb.8-r8.3`.

Report issues [on GitHub](https://github.com/smash-transport/smash/issues) or contact us by  [✉️ email](mailto:elfner@itp.uni-frankfurt.de).

## How to build and install SMASH

In the following you can find a minimal quick start guide.
Refer to the [INSTALL](INSTALL.md) file for more detailed information.

### Prerequisites

SMASH is known to compile and work on 64-bit little endian machines (most CPUs are such) with UNIX-like operating systems (e.g. GNU/Linux, MacOS) and one of the following compilers (which have the required C++17 features):

| Compiler   | Required version |
|  :---:     |       :---:      |
| GCC        |  8.0 or higher   |
| Clang      |  7.0 or higher   |
| Apple clang| 11.0 or higher   |

Any different operating system and/or compiler and/or endianness is not officially supported and SMASH will ask you to continue at your own risk before compilation.

SMASH requires the following tools and libraries:

| Software | Required version |
|  :---:   |       :---:      |
| [CMake](https://cmake.org) | 3.16 or higher |
| [GNU Scientific Library (GSL)](https://www.gnu.org/software/gsl/) | 2.0  or higher |
| [Eigen3 library](http://eigen.tuxfamily.org) | 3.0  or higher |
| [Pythia](https://pythia.org) | 8.316 |

Support for ROOT, HepMC3 and Rivet output is automatically enabled if a suitable version is found on the system:

| Software | Required version |
|  :---:   |       :---:      |
| ROOT     | 5.34 or higher   |
| HepMC3   | 3.2.3 or higher  |
| Rivet    | 3.1.4 or higher  |

### Compilation and installation

Clone the SMASH repository with the help of `git clone`.
From within the cloned repository, use the following commands to build the codebase in a separate directory:
```console
mkdir build
cd build
cmake -DPythia_CONFIG_EXECUTABLE=/path/to/pythia8316/bin/pythia8-config ..
make
```
Please note that the `make` command builds everything (executables, tests, and libraries) and this might take a while.
You can use `make smash` if you are interested only in the SMASH executable or use `make smash_shared` to exclusively build the libraries (needed e.g. in another project using SMASH as a library -- refer to the [INSTALL](INSTALL.md) FAQ for more detailed information).

You can run SMASH with specific settings (e.g. at a given collision energy or impact parameter) by modifying the config.yaml file, for example with
```console
vi config.yaml
./smash
```
Refer to the [section below](#running-smash-with-example-input-files) for more information.

If you want to install SMASH system-wide (into `/usr/local`) use
```console
make install
```

⚠️ **NOTE:** All commands above are the bare minimum needed for an installation.
It is not guaranteed that this minimum setup is appropriate for your needs or your specific computing environment.
For example, several different options can be passed e.g. to the `cmake` command.
We strongly advise you to further refer to the [INSTALL](INSTALL.md) file for more guidance, especially if you encounter any issues.


## Using the Docker containers

As an alternative to building or installing SMASH, a Docker image of the latest or recently tagged version can be pulled from the Github container registry.
Get the newest version with
```console
docker pull ghcr.io/smash-transport/smash:newest
```

Start the container with
```console
docker run -it ghcr.io/smash-transport/smash:newest
```

A ready-to-use executable of SMASH is found in the `smash_bin` directory.
Run it as explained below.
If needed, SMASH can also be built inside the container as explained in the previous section (the SMASH source files and Pythia are also found in the `/SMASH` directory).

Two container versions of SMASH are offered:
* a small version (`ghcr.io/smash-transport/smash`) with a minimal set of dependencies
pre-installed and
* a large version with all possible external dependencies, e.g. ROOT, HepMC and Rivet, already included (`ghcr.io/smash-transport/smash-max`).

Note that running SMASH inside of a Docker container might negatively affect performance.
More information about containers usage can be found [here](containers/README.md).

#### Note for users with ARM CPUs (e.g. Apple M1/M2 chips)

Our Docker images are prepared for the x86-64 CPU architecture.
To make them compatible with computers with ARM CPUs (like in the case of Apple M1 and M2 chips),
`docker` must be launched with the `--platform=linux/amd64` option.
For example:
```console
docker run --platform=linux/amd64 -it ghcr.io/smash-transport/smash:newest
```
However, this is not always guaranteed to work and it might be necessary to build an image for the ARM architecture, as described [here](containers/README.md).

<a id="running-smash-with-example-input-files"></a>

## Running SMASH with Example Input Files

SMASH ships example configuration files for running in the collider, box, sphere, and list mode (`Modus` in the configuration jargon).
By default, i.e. by running `./smash`, the simulation is set up from the collider configuration file, called _config.yaml_, and using the default particles and decay modes files (_particles.txt_ and _decaymodes.txt_, respectively).
They are located in the repository ***input*** folder.

Additionally, example configuration files for the box, sphere and list modus can be found in the respective directories ***input/{box,sphere,list}***.
If needed, e.g. in the case of a box simulation, different default particles and decay modes files can be used.
Examples for these are also provided in ***input/box***.

Finally, for the list modus, an input list file to be read in is required and an example is provided as _input/list/example_list0_.

In general, to run SMASH with a non-default configuration file, use the `-i` command.
For example, for the sphere or list example file, from the ***build*** folder, use:
```console
    ./smash -i ../input/sphere/config.yaml
    ./smash -i ../input/list/config.yaml
```

Furthermore, if using non-default particles and decay modes files is necessary, these can be specified through the `-p` and `-d` options.
In the box and the dileptons example, always from the ***build*** folder, this means:
```console
./smash -i ../input/box/config.yaml -p ../input/box/particles.txt -d ../input/box/decaymodes.txt
./smash -i ../input/dileptons/config.yaml -d ../input/dileptons/decaymodes.txt
```

All available command line options for SMASH can be viewed with
```console
./smash -h
```
To run SMASH completely silently for production runs, we recommend to suppress the standard output via e.g.
```console
./smash > /dev/null
```
and it might be useful to redirect warnings and error messages, that will still be displayed, to a file:
```console
./smash > /dev/null 2> /path/to/error-and-warnings-file
```


## License

SMASH is licensed under the terms of the GNU General Public License, Version 3 or above.
The build scripts are licensed under terms of the BSD 3-clause license.
For more information, see [LICENSE](LICENSE).


## Projects Using SMASH

SMASH source and documentation are provided to check and reproduce published results of the authors.
Cooperation and joint projects with outside researchers are encouraged and comparison to results by experimental collaborations is supported.
If you are interested in starting a project, please contact us to avoid interference with current thesis topics.
If your project involves changes to the code, please refer to [CONTRIBUTING](CONTRIBUTING.md) for coding guidelines and helpful tools.
SMASH can also be used as a 3rd party library -- refer to the [INSTALL](INSTALL.md) FAQ for more detailed information.
