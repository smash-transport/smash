# syntax=docker/dockerfile:1.4

FROM ubuntu:22.04

# Set bash as the default shell in this Dockerfile

SHELL ["/bin/bash", "-c"]

# Set the proper general cpu family architecture flag.
# Please, for more information refer to the documention of your compiler.

ARG TARGET_ARCHITECTURE="x86-64"

# Add some info about authors and contents

LABEL maintainer.name="SMASH team"
LABEL maintainer.website="https://smash-transport.github.io"
LABEL description="Basic SMASH image (sources, compiler and bare executable, without ROOT, HepMC3 and Rivet)."
LABEL org.opencontainers.image.source https://github.com/smash-transport/smash

# Install minimal set of requirements

RUN apt-get update && \
apt-get -y upgrade && \
DEBIAN_FRONTEND=noninteractive apt-get -y install \
bash \
cmake \
git \
g++ \
libeigen3-dev \
libgsl-dev \
rsync \
wget \
vim

# Set cpu architecture flag for compilation

ARG C_CXX_ARCH_FLAG="-march=${TARGET_ARCHITECTURE}"

# Create main project directory

WORKDIR /SMASH

# Create directory for mounting filesystems

RUN mkdir /MNT

# Install pythia

ARG pythiaV="8310"
ENV PYTHIA8=/SMASH/pythia_install/
RUN wget https://pythia.org/download/pythia83/pythia${pythiaV}.tgz && \
tar xf pythia${pythiaV}.tgz && rm pythia${pythiaV}.tgz && \
cd pythia${pythiaV} && \
./configure --cxx-common="${C_CXX_ARCH_FLAG} -O3 -std=c++17 -fPIC -pthread" --prefix=${PYTHIA8} && \
make -j$(nproc) && \
make install && \
echo ${PYTHIA8}/lib >> /etc/ld.so.conf && ldconfig

ENV PYTHIA8DATA=$PYTHIA8/share/Pythia8/xmldoc
ENV PYTHIA_LIBRARY=$PYTHIA8/lib
ENV PYTHIA_INCLUDE=$PYTHIA8/include
ENV LD_LIBRARY_PATH=$PYTHIA_LIBRARY:$LD_LIBRARY_PATH
ENV PATH=$PYTHIA8/bin:$PATH

# Get and build public SMASH

RUN cd /SMASH && git clone https://github.com/smash-transport/smash.git && \
mkdir -p smash/build && \
mkdir -p smash_bin && \
cd smash/build && \
cmake .. -DPythia_CONFIG_EXECUTABLE=/SMASH/pythia${pythiaV}/bin/pythia8-config -DCMAKE_C_FLAGS="${C_CXX_ARCH_FLAG}" \
-DCMAKE_CXX_FLAGS="${C_CXX_ARCH_FLAG}" || echo "cmake SMASH build configuration failed!"
RUN cd /SMASH/smash/build && make smash -j$(nproc) || echo "make SMASH compilation failed!"
RUN <<EOF
if [ -f /SMASH/smash/build/smash ]
then
  cp /SMASH/smash/build/smash /SMASH/smash_bin/
  cp /SMASH/smash/build/config.yaml /SMASH/smash_bin/
  rm -r /SMASH/smash/build
else
  rm -r /SMASH/smash /SMASH/smash_bin
fi
EOF

# Set PATH

ENV PATH=/SMASH/smash_bin:$PATH
