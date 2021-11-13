# README From Docker to Singularity

In most cases Singularity is able to import or directly run commands from Docker containers without significant problems.
In both situations Singularity caches information in a subdirectory of ~/.singularity.
More information at: [https://sylabs.io/guides/3.6/user-guide/singularity_and_docker.html](https://sylabs.io/guides/3.6/user-guide/singularity_and_docker.html).

### Import a Docker container from an online registry

It is possible to retrieve a Docker container _docker-cont_ from an online registry _ghcr.io/repositoryname/_ and tranform it into a Singularity container, named in this case _mycontainer.sif_, with:

```
singularity pull mycontainer.sif docker://ghcr.io/repositoryname/docker-cont
```

It is possible also to execute single commands, for example:

```
singularity exec docker://ghcr.io/repositoryname/docker-cont bash
```
launches a bash shell within the container. If the Docker image has been already cached, Singularity uses the local copy in .sif format, otherwise it downloads the Docker image and converts it into .sif format in the internal cache without creating another additional .sif file. Then, of course, Singularity executes the command.


### Import a local Docker container

First, the local Docker container should be exported as a tar archive. For example, assuming that we want to export the image _abcd1234_:

```
docker save --output=abcd1234-docker-cont.tar abcd1234
```

In general it is convenient to compress the tar archive before transmission over internet and decompress it before conversion.

To convert _abcd1234-docker-cont.tar_ into _mycontainer.sif_ (both in the same current working directory):

```
singularity build mycontainer.sif docker-archive://abcd1234.tar
```

This kind of operation does not require root privileges.
