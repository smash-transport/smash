# README Singularity

#### Contents of this README
- Quick build instructions
- Short introduction about containers and Singularity
- More extended build instructions and explanations

If you are searching information about how to build a container with Docker, please, give a look at [README_Docker.md](README_Docker.md).

## Quick build instructions

The following procedure has been tested on an Ubuntu 20.04 GNU/Linux distribution with a Bash shell and Singularity version 3.6.4, but it should work without any change also in other environments running a similar version of Singularity.

First of all, we assume that you:
- have installed Singularity (if not, there are a few hints later in this document) and you can invoke it from the command line of a shell
- you have _root_ or _sudo_ user access on the machine used to build the container
- are connected to internet to download additional software/packages

In the source directory containing the files: _smash_v2_0__ubuntu_20_04.def_ and _list_of_packages_for_ubuntu20_04.txt_ issue the command:

`sudo singularity build smash_v2_0__ubuntu_20_04.sif smash_v2_0__ubuntu_20_04.def`

If everything goes well, after some time (more than one hour on an i7-2670QM with a 16Mb/s internet connection) you will get _smash_v2_0__ubuntu_20_04.sif_.

## Short introduction about containers and Singularity
The purpose of containers is to provide a complete environment, with all what is needed to run a certain software, with a high degree of independence from the system in which the container is used, so that the user does not have to worry about the native version of the OS, the libraries and so on of the hosting system, making life easier.

The most common container system is [Docker](https://www.docker.com/), but it has some characteristics that don't make it suitable for High Performance Computing systems (HPC), like the management of the filesystem or the user privileges. [Singularity](https://sylabs.io/) has been developed to address these issues, providing a flexible, secure and highly performant framework for using containers in HPC Clusters.

To get more information about what are containers and Singularity, we suggest to watch these videos on youtube:
* [https://www.youtube.com/watch?v=DA87Ba2dpNM](https://www.youtube.com/watch?v=DA87Ba2dpNM), a slightly old, but nice introduction by Gregory M. Kurtzer, providing an overall good picture of why containers are useful and what are the needs that Singularity tries to fullfil and that other systems, like the more famous Docker, are not able to satisfy
* [https://www.youtube.com/watch?v=ERxK0i0SDA0](https://www.youtube.com/watch?v=ERxK0i0SDA0), a more practical hands on about building a container

This is a list of written references:
* [https://sylabs.io/guides/3.6/user-guide/](https://sylabs.io/guides/3.6/user-guide/), the official user guide of Singularity
* [https://hpc.gsi.de/virgo/containers/overview.html](https://hpc.gsi.de/virgo/containers/overview.html), the guide about Singularity and containers in the official documentation of the Virgo Cluster

At the moment of writing this notes, there is an active development of Singularity and features and standard/conventions can change from one version to the other. We recommend to always check the official documentation of the version that you are using and to build the container with the same version of Singularity running on the Virgo Cluster (or of the machine in which you will deploy it).

## More extended build instructions and explanations

### How to build the Singularity container for SMASH
Most of the times containers under Singularity are run as normal users and actually one of the most interesting security features of Singularity is that it prevents escalation of privileges. However, building a container requires root privileges and for normal users it has to be done outside of a Cluster.

Here are the instructions about how to install Singularity on Linux: [https://sylabs.io/guides/3.6/admin-guide/installation.html](https://sylabs.io/guides/3.6/admin-guide/installation.html). The compilation requires a few dependences and [Go](https://golang.org/), but the instructions are clear and it should be pretty straightforward to obtain a working Singularity installation. If Singularity is not installed in a directory included in the PATH environment variable, i.e. in the list of directories in which GNU/Linux searches the executables, it might be convenient to add it, e.g.: `export PATH=/opt/singularity/bin:$PATH`.

### Building a container from a definition file
It is possible to build a container, let's call it _container.sif_, by writing the procedure in a file, let's call it _definitions.def_, and then issue, as __root__ user (or with sudo):

`singularity build container.sif definitions.def`

To get some basic informations about a container called _container.sif_:
* `singularity inspect container.sif`
* `singularity run-help container.sif`

It is possible to extract the definitions used to build a container with the command:

`singularity sif dump 1 container.sif >> definitions.def`

The ">>" redirects the output to the file _definitions.def_.

The file with the definitions is composed by several sections (for example %post, %environment, %startscript, and so on) that instruct Singularity not only about how to build the container, but also how to set up the environment when it is launched, what to execute when it starts, which help message to print, and so on. We refer to [https://sylabs.io/guides/3.6/user-guide/definition_files.html](https://sylabs.io/guides/3.6/user-guide/definition_files.html) for more details.

Here we briefly summarize the steps taken to build the SMASH container, that one should be able to identify in the definition file (_.def_  extension) with the instructions.
The def file is divided into several section, delimited by labels like %setup, %files, %post, %environment, %label and %help. We refer to the official Singularity documentation for further details.
* we start by reusing a small container which provides a basic GNU/Linux installation.
For example, one can directly download from the official library repository of Singularity a basic version of Debian 10:
```
Bootstrap: library
From: debian
```
However, when using a Debian 10 the final container showed some memory leak problems with root when testing SMASH, therefore, for the first official container, we used an Ubuntu 20.04 distribution, taken from the official Docker repository with:
```
Bootstrap: docker
From: ubuntu:20.04
```
It is also possible to download an image with:
`singularity pull docker://ubuntu`
and then use the local file (which has been automatically converted in sif format by Singularity when retrieved) with:
```
Bootstrap: localimage
From: smash_stuff/ubuntu_latest.sif
```
* the %files section, often not used when preparing containers for SMASH, could be exploited to copy in the container some local files. Keeping a local copy of some software might be a good idea to save time in case of multiple container builidings, avoiding to download many times the same file, or to be sure that a certain version, working perfectly, will always be available also in the future.
* we choose some local (e.g. located in Germany) repositories of the distribution packages, we update the base system and we install several other packages
* in the %post section we download and compile SMASH and several auxiliary libraries/programs
* we modify one of the final startup scripts inside the _/.singularity.d/_ directory that will be executed when the container starts to set up the prompt (this might be a temporary workaround, but in principle this approach might be used also for other purposes, albeit the preferred method to setup the environment is to use the _%environment_ section and let Singularity to take care to prepare the startup scripts)
* in the _%environment_ section we define several environment variables that will be used during the execution of the container
* we fill in some basic information in the _%label_ and _%help_ sections

Using this procedure, all the needed files are downloaded when building the container. However, this will be done every time that one builds the container, therefore, in some cases, it might be more convenient to collect some stuff in local directories and then copy their content into the container by including a _%files_ section in the _definitions.def_ (as the GSI IT Administrators have done for the official Virgo containers), with a list of sources in the local filesystem and destinations in the future filesystem of the container (see the Singularity documentation for more details).

### Building a container using a sandbox

It is possible to unpack a container into a directory, modify it interactively and then pack again the results into a new container. For example, one can:
* create a temporary directory: `mkdir container_dir_sandbox`
* unpack the container _container.sif_ inside the directory: `singularity build -s container_dir_sandbox container.sif`
* run a shell inside the container in writing mode: `singularity shell -w container_dir_sandbox`
* use the shell to modify the software inside the container
* build an updated version of the container (_new_container.sif_) with: `singularity build new_container.sif container_dir_sandbox`

### Basic usage of the container:

You can run a shell within the container with: `singularity shell container.sif` where _container.sif_ is actually the filename of the container that you want to use (e.g. _smash_v2_0__ubuntu_20_04.sif_).
This shell sees as system directories the directories created inside the container, which are all read-only, but the directories that are not part of the container, like the home directory of the user, are still accessible with the same permissions granted to the user by the host system.
The _smash_ executable is located in the /SMASH/smash_bin directory, while /SMASH/smash contains the sources.
If the user wants to run the tests, she/he should create a subdirectory in a directory in which she/he has writing permissions `mkdir -p bb && cd bb` compile SMASH there `cmake /SMASH/smash -DCMAKE_INSTALL_PREFIX=/SMASH/hepmc3-install/`, `make -j4` and then run the tests `ctest -j4`. The number after the -j option tells how many threads you want to use, very useful if you have a cpu with multiple cores or, at least, hyperthreading.

### Final comments:

We recall that, apart from SMASH, there are many official containers ready to use. For example, one can download a CentOS distribution (e.g. to test if the compilation of a program succeed or fails with a certain version of gcc/libraries) from the Singularity library with:

`singularity pull --arch amd64 library://library/default/centos:8`

We also recall that there is the chance to use/reuse Docker containers, too. Please, have a look at: [https://sylabs.io/guides/3.6/user-guide/singularity_and_docker.html](https://sylabs.io/guides/3.6/user-guide/singularity_and_docker.html) for further information. Of course, one can also install Docker and take directly advantage of the many ready-to-use containers for this framework. Almost always Docker containers can be imported in Singularity, the pull command will convert them into the native _sif_ format. So, it is also possible to build a container for SMASH with Docker and then convert its format into that used by Singularity. If are interested in this approach, please, give a look to [README_Docker.md](README_Docker.md).

Please, note that this guide does not mention many important features and technical details that you can find in the official user and admin guides of Singularity.
