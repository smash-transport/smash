FROM ubuntu:20.04

# Install minimal set of requirements
RUN apt-get update && \
apt-get -y upgrade && \
DEBIAN_FRONTEND=noninteractive apt-get -y install \
cmake \
git \
g++ \
libboost1.71-dev \
libboost-filesystem1.71-dev \
libeigen3-dev \
libgsl-dev \
rsync \
wget \
vim

# Create main project directory
WORKDIR /SMASH

# Install pythia
ARG pythiaV="8307"
RUN wget https://pythia.org/download/pythia83/pythia${pythiaV}.tgz && \
tar xf pythia${pythiaV}.tgz && rm pythia${pythiaV}.tgz && \
cd pythia${pythiaV} && \
./configure --cxx-common='-std=c++11 -mfpmath=sse -O3 -fPIC -march=x86-64' && \
make -j$(nproc) && \
make install && \
echo /SMASH/pythia${pythiaV}/lib >> /etc/ld.so.conf && ldconfig

# Get and build public SMASH
RUN git clone https://github.com/smash-transport/smash.git && \
mkdir -p smash/build && \
mkdir -p smash_bin && \
cd smash/build && \
cmake .. -DPythia_CONFIG_EXECUTABLE=/SMASH/pythia${pythiaV}/bin/pythia8-config && \
make smash -j$(nproc) && \
cp smash  /SMASH/smash_bin/ && \
cp config.yaml /SMASH/smash_bin/ && \
rm -r /SMASH/smash/build
