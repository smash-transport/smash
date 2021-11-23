# README Docker

Information about Docker can be found at: https://docs.docker.com/. The links to download Docker for Linux, Mac and Windows10 can be found at: https://www.docker.com/get-started.
The instructions to build SMASH are written in the file `Dockerfile` and, apart from the different syntax, the procedure is the same as in the case of Singularity (see [README_Singularity.md](README_Singularity.md)). As explained in the [README.md](../README.md), Docker images are provided as part of the Github repository and can be pulled from there and run directly

### Building a Docker image for SMASH

Assuming that Docker is already installed, one can also build images locally by executing in a terminal, in the same directory of `Dockerfile`:

```
docker build .
```

Docker keeps the information about the container in a directory, e.g., under GNU/Linux, in /var/lib/docker/, instead of just in a file like Singularity. One can get an overview of the current Docker containers/installation with: `docker info` and, in particular, it is possible to list the available images with: `docker images`
The images have each its own id. One can recognize the last one just produced by the creation date (tags are also possible, but not discussed here). Run a container with `docker run -it <image_id>`. You might need to use `sudo` in front of the docker commands.
To export an image into a file: `docker save -o name_of_Docker_container_image.tar image_id`. It is suggested to compress the file (e.g. with gzip) before uploading it into the cluster or somewhere else.
Once the tar archive of the Docker container is in the cluster, one can transform it into a singularity container with: `singularity build name_of_the_new_singularity_container.sif docker-archive://name_of_Docker_container_image.tar`

We recall that, to change the default prompt of the bash shell, the user can set up the environment variable _PS1_ (see [https://www.gnu.org/software/bash/manual/bash.html#Controlling-the-Prompt](https://www.gnu.org/software/bash/manual/bash.html#Controlling-the-Prompt) for more information).

Note that default resources allocated for the Docker deamon might not be enough. In particular, the virtual memory might be insufficent and lead to crashes. Try increasing the resoruces, if the build process gets terminated.

### Using a Docker container for development

If one wants to use the container environment for development, it is useful to mirror a smash development repository into a container directory with running the container with `-v` option as follows:

```
docker run -it -v path/to/smash/repo:/SMASH/smash_local  <image_id or tag>

```
This creates the directory `/SMASH/smash_local` which matches the smash directory. All local changes will be reflected in this directory and SMASH can be build with those changes in the container.
