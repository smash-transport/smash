## Docker containers

Information about Docker can be found on the [Docker official page](https://docs.docker.com/).
The links to download Docker for Linux, Mac and Windows10 are also [officially available](https://www.docker.com/get-started).
The instructions to build SMASH are written in the file _Dockerfile_ and, apart from the different syntax, the procedure is the same as in the case of Singularity (see [following section](#docker-to-singularity)).
As explained in the [SMASH README file](../README.md), Docker images are provided as [packages in the Github organisation](https://github.com/orgs/smash-transport/packages) and can be pulled from there and directly run. 

### Building a Docker image for SMASH

Assuming that Docker is already installed, one can also build images locally by executing
```console
docker buildx build -f Dockerfile .
```
in a terminal, in the same directory of _Dockerfile_.
Our images offer the possibility to choose the architecture to target at compilation time in the following way:
```console
docker buildx build -f Dockerfile --build-arg="TARGET_ARCHITECTURE=x86-64" .
```
Of course, `x86-64` can be changed using any general cpu family architecture flag.

The option `-t <tag_name>` can be added to assign the `<tag_name>` tag to the built image.

#### Note for users with ARM CPUs (e.g. Apple M1/M2 chips)

Our Docker images are by default prepared for the x86-64 CPU architecture.
To make them compatible with computers with ARM CPUs (like in the case of Apple M1 and M2 chips),
`docker` must be launched with the `--platform=linux/amd64` option.
For example
```console
docker buildx build -f Dockerfile .
```
becomes
```console
docker buildx build --platform=linux/amd64 -f Dockerfile .
```
However, in this way, roughly speaking, Docker adds a layer of virtualization and an image for an x86-64 machine is built (and potentially run with the same option) e.g. from an arm64 one.
This can lead to problems and it is not guaranteed to work, though.

If you plan to use the image on ARM CPUs **only**, then you do not need to add any further virtualization layer.
Instead, you can use the `--build-arg` option of the `docker` command to target your architecture at compilation time and build the image from and for your ARM CPU.
For instance
```console
docker buildx build -f Dockerfile --build-arg="TARGET_ARCHITECTURE=native" .
```
will optimize the image for the architecture of the machine on which the image itself is built.

**NOTE:** To build the "max" SMASH image on your ARM machine, you need to first build the "min" one and specify it in the `FROM` command at the beginning of the "max" image Dockerfile, using there the tag name used to tag the "min" image.

#### Docker in a nutshell

Docker keeps the information about the container in a directory (e.g. in ***/var/lib/docker/*** under GNU/Linux), instead of just in a file like Singularity.
One can get an overview of the current Docker containers/installation with `docker info` and, in particular, it is possible to list the available images with `docker images`.
The images have each its own id.
One can recognize the last one just produced by the creation date (tags are also possible, but not discussed here).
Run a container with `docker run -it <image_id>`.
You might need to use `sudo` in front of the docker commands.
To export an image into a file, use `docker save -o name_of_Docker_container_image.tar image_id`.
It is suggested to compress the file (e.g. with gzip) before uploading it into the cluster or somewhere else.
Once the tar archive of the Docker container is in the cluster, one can transform it into a singularity container via
```console
singularity build name_of_the_new_singularity_container.sif docker-archive://name_of_Docker_container_image.tar
```

We recall that, to change the default prompt of the bash shell, the user can set up the environment variable `PS1` (refer to [the bash manual](https://www.gnu.org/software/bash/manual/bash.html#Controlling-the-Prompt) for more information).

Note that default resources allocated for the Docker deamon might not be enough. In particular, the virtual memory might be insufficent and lead to crashes.
Try increasing those resoruces, if the build process gets suddenly terminated.

#### Using a Docker container for development

If you want to use the container environment for development, it is useful to mirror a smash development repository into a container directory with running the container with `-v` option as follows:
```console
docker run -it -v path/to/smash/repo:/SMASH/smash_local  <image_id or tag>
```
This creates the directory `/SMASH/smash_local` which matches the smash directory.
All local changes will be reflected in this directory and SMASH can be build with those changes in the container.
However, sometimes writing error permissions might occur, usually because the container user ID differs from the host one.
You can make them match via
```console
docker run -it -u $(id -u) -v path/to/smash/repo:/SMASH/smash_local <image_id or tag>
```
and you can e.g. refer to [this SO answer](https://stackoverflow.com/a/66350210/14967071) for more information.



<a id="docker-to-singularity"></a>

## From Docker to Singularity/Apptainer

We remind that in many HPC systems Singularity is actually replaced by Apptainer, but the two projects have tight contacts and an extremely high compatibility.
In the examples we will refer only to Singularity, but the commands are exactly the same also with Apptainer, including the name of the executable `singularity`.
In most cases Singularity and Apptainer are able to import or directly run commands from Docker containers without significant problems.
In both situations Singularity and Apptainer cache information in a subdirectory of ***~/.singularity*** or ***~/.apptainer***, respectively.
More information can be found [here](https://docs.sylabs.io/guides/latest/user-guide/singularity_and_docker.html) and [here](https://apptainer.org/docs/user/latest/docker_and_oci.html).

### Import a Docker image from an online registry

It is possible to retrieve a Docker image _docker-image_ from an online registry `ghcr.io/repositoryname/` and tranform it into a Singularity image, named in this case _my-singularity-image.sif_, with:
```console
singularity pull my-singularity-image.sif docker://ghcr.io/repositoryname/docker-image
```

For example, in the case of the official SMASH basic image:
```console
singularity pull smash.sif docker://ghcr.io/smash-transport/smash:newest
```

It is possible also to execute single commands in the container.
For example
```console
singularity exec docker://ghcr.io/repositoryname/docker-image bash
```
launches a bash shell within the container.
If the Docker image has been already cached, Singularity uses the local copy in `.sif` format, otherwise it downloads the Docker image and converts it into `.sif` format in the internal cache without creating another additional `.sif` file.
Then, of course, Singularity executes the command.

### Import a local Docker image

First, the local Docker image should be exported as a tar archive.
For example, assuming that we want to export the image _abcd1234_, use
```console
docker save --output=abcd1234-docker-cont.tar abcd1234
```

In general it is convenient to compress the tar archive before transmission over internet and decompress it before conversion.

To convert _abcd1234-docker-cont.tar_ into _mycontainer-image.sif_ (both in the same current working directory) use
```console
singularity build mycontainer-image.sif docker-archive://abcd1234.tar
```

This kind of operation does not require root privileges.
