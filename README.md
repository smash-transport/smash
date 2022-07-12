# SMASH

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3484711.svg)](https://doi.org/10.5281/zenodo.3484711)

SMASH (Simulating Many Accelerated Strongly-interacting Hadrons) is a relativistic hadronic transport approach for the dynamical description of heavy-ion reactions.
Please see [Phys. Rev. C 94, 054905 (2016)](https://arxiv.org/abs/1606.06642) for details and cite this reference together with the [software DOI](https://doi.org/10.5281/zenodo.3484711) for the specific code version employed, if you are using SMASH.
A BibTeX entry for the software DOI is found on the respective Zenodo pages.

See [CONTRIBUTING](CONTRIBUTING.md) for development hints.
A complete [User Guide](https://theory.gsi.de/~smash/userguide/current/) as well as a more detailed [development documentation](http://theory.gsi.de/~smash/doc/current/) are available.

If Pythia is used, please cite the following references: T. Sjöstrand, S. Mrenna and P. Skands, [JHEP05 (2006) 026](https://arxiv.org/abs/hep-ph/0603175) and [Comput. Phys. Comm. 178 (2008)](https://arxiv.org/abs/0710.3820).

Report issues [on GitHub](https://github.com/smash-transport/smash/issues) or contact us by email at elfner@itp.uni-frankfurt.de.

## How to build and install SMASH

In the following you can find a minimal quick start guide.
Refer to the [INSTALL](INSTALL.md) file for more detailed information.

### Prerequisites

SMASH is known to compile and work on little endian machines with UNIX-like operating systems (e.g. GNU/Linux, MacOS) and one of the following compilers (which have the required C++17 features).

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
| [Pythia](https://pythia.org) | 8.307 |

Support for ROOT, HepMC3 and Rivet output is automatically enabled if a suitable version is found on the system.

| Software | Required version |
|  :---:   |       :---:      |
| ROOT     | 5.34 or higher   |
| HepMC3   | 3.2.3 or higher  |
| Rivet    | 3.1.4            |

### Compilation and installation

From within the SMASH repository, use the following commands to build the codebase in a separate directory:
```console
mkdir build
cd build
cmake -DPythia_CONFIG_EXECUTABLE=/path/to/pythia8307/bin/pythia8-config ..
make
```

You can run SMASH with specific settings (e.g. at a given collision energy or impact parameter) by modifying the config.yaml file, for example with
```console
vi config.yaml
./smash
```
Refer to the [section below](README.md#running-smash-with-example-input-files) for more information.

If you want to install SMASH system-wide (into `/usr/local`) use
```console
make install
```

**NOTE:** All commands above are the bare minimum needed for an installation.
It is not guaranteed that this minimum setup is appropriate for your needs or your specific computing environment.
For example, several different options can be passed e.g. to the `cmake` command.
We strongly advise you to further refer to the [INSTALL](INSTALL.md) file for more guidance, especially if you encounter any issues.


## Using the Docker containers

Alternatively to building or installing SMASH, a Docker image of the latest or recently tagged version can be pulled from the Github container registry.
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
If needed, SMASH can also be build inside the container as explained in the previous section (the SMASH source files and Pythia are also found in the `/SMASH` directory).

Two container versions of SMASH are offered: a small version (`ghcr.io/smash-transport/smash`) with a minimal set of dependencies
pre-installed and a maximum version with all possible external dependencies, e.g. ROOT, HepMC and Rivet, already included (`ghcr.io/smash-transport/smash-max`).
Running SMASH inside of a Docker container might negatively affect performance.
More information on container usage is found in the README files in the `containers` directory.


## Running SMASH with Example Input Files

SMASH ships with example configuration files for the collider, box, sphere, and list modus.
By default, i.e. by running `./smash`, the simulation is set up from the collider configuration file, called `config.yaml`, and using the default particles and decay modes files (`particles.txt` and `decaymodes.txt`, respectively).
They are located in the repository `input` folder.

Additionally, example configuration files for the box, sphere and list modus can be found in the respective directories `input/{box,sphere,list}`.
In case of the box simulation, in order to allow for equilibration, different default particles and decay modes files need to be used.
These files are also provided in `input/box`.
For the list modus, an input list file to be read in is required.
This file, `example_list0`, is located in `input/list`.

To run SMASH with a non-default configuration file, use the `-i` command.
For example, for the sphere or list example file, from the `build` folder, use:
```console
    ./smash -i ../input/sphere/config.yaml
    ./smash -i ../input/list/config.yaml
```

Furthermore, using non-default particles and decay modes files is necessary, and these can be specified through the `-p` and `-d` options.
In the box and the dileptons example, always from the `build` folder, this means:
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

Warnings and error messages will still be displayed.


## License

SMASH is licensed under the terms of the GNU General Public License, Version 3 or above.
The build scripts are licensed under terms of the BSD 3-clause license.
For more information, see [LICENSE](LICENSE).


## Projects Using SMASH

SMASH source and documentation are provided to check and reproduce published results of the authors.
Cooperation and joint projects with outside researchers are encouraged and comparison to results by experimental collaborations is supported.
SMASH can be used as a 3rd party library, for examples see the `examples` folder in the repository.
If you are interested in starting a project, please contact us to avoid interference with current thesis topics.
If your project involves changes to the code, please refer to [CONTRIBUTING](CONTRIBUTING.md) for coding guidelines and helpful tools.
