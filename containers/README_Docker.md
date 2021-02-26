# README Docker

Information about Docker can be found at: https://docs.docker.com/. The links to download Docker for Linux, Mac and Windows10 can be found at: https://www.docker.com/get-started.
The instructions to build SMASH are written in the file `Dockerfile` and, apart from the different syntax, the procedure is the same as in the case of Singularity (see [README_Singularity.md](README_Singularity.md)):
- download a basic GNU/Linux distribution from a Docker repository
- choose some distribution mirros, update the distribution and add several additional packages
- use pip to install python2, not officialy supported anymore
- download, compile and install several programs and libraries used by SMASH
- download and compile SMASH
- export some environment variables that will be visible when starting/using the container

The number of processes in concurrent compilation is set in the line:

```
ARG NUMBER_OF_PROCS_FOR_COMPILATION=2
```

located after the list of packages to upgrade the basic GNU/Linux distribution. This value can be increased if the user has more than 2 cores available. Please, beware that under Mac and Windows Docker exploits a virtual machine that, by default, has only half of the physical cpu cores and 2GB of ram, therefore it is necessary to modify accordingly also these settings (please, refer to the official Docker documentation to see how).

### Building a Docker container for SMASH

Assuming that Docker is already installed, one has to execute in a terminal, in the same directory of `Dockerfile`:

```
sudo docker build -f Dockerfile .
```

Docker keeps the information about the container in a directory, e.g., under GNU/Linux, in /var/lib/docker/, instead of just in a file like Singularity.
One can get an overview of the current Docker containers/installation with: `sudo docker info` and, in particular, it is possible to list the available images with: `sudo docker images`
The images have each its own id. One can recognize the last one just produced by the creation date (tags are also possible, but not discussed here).
To export an image into a file: `sudo docker save -o name_of_Docker_container_image.tar image_id`. It is suggested to compress the file (e.g. with gzip) before uploading it into the cluster or somewhere else.
Once the tar archive of the Docker container is in the cluster, one can transform it into a singularity container with: `singularity build name_of_the_new_singularity_container.sif docker-archive://name_of_Docker_container_image.tar`

The creation of a container with Docker has been tried on Linux (Ubuntu 20.04) and with Windows10 Home ed. (with WSL2 installed), in the second case in a PowerShell without using _sudo_.

We recall that, to change the default prompt of the bash shell, the user can set up the environment variable _PS1_ (see [https://www.gnu.org/software/bash/manual/bash.html#Controlling-the-Prompt](https://www.gnu.org/software/bash/manual/bash.html#Controlling-the-Prompt) for more information).
