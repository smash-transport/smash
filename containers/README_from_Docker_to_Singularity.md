# README From Docker to Singularity

In most cases Singularity is able to import and use Docker containers without significant problems.
Singularity can also directly run commands in Docker containers, but, since a cached converted copy is stored in a subdirectory of ~/.singularity anyway, in many circumstances it may be better to explicitly perform a conversion into the native Singularity .sif format.
More information at: [https://sylabs.io/guides/3.6/user-guide/singularity_and_docker.html](https://sylabs.io/guides/3.6/user-guide/singularity_and_docker.html).

### Import a Docker container from an online registry

It is possible to retrieve a Docker container _docker-cont_ from an online registry _ghcr.io/developername/_ and tranform it into a Singularity container, named in this case _mycontainer.sif_, with:

```
singularity pull mycontainer.sif docker://ghcr.io/developername/docker-cont
```

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
