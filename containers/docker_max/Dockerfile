# syntax=docker/dockerfile:1.4

# Import the base smash Docker container image.
# Adjust the address as needed. In particular, be sure to use a base image compatible with your cpu.
# The official Docker images are prepared for a generic Intel/AMD x86-64 cpu.

FROM ghcr.io/smash-transport/smash:newest

# Set bash as the default shell in this Dockerfile

SHELL ["/bin/bash", "-c"]

# Set the proper general cpu family architecture flag.
# Please, for more information refer to the documentation of your compiler.

ARG TARGET_ARCHITECTURE="x86-64"

# Declare all together the versions of the going to be installed codebases

ARG ROOT_VERSION="6.34.00"
ARG HEPMC_VERSION="3.2.6"
ARG CPPCHECK_VERSION="2.8"
ARG CPPLINT_VERSION="1.6.0"
ARG RIVET_VERSION="3.1.7"
ARG YODA_VERSION="1.9.7"
ARG FASTJET_VERSION="3.4.0"
ARG FJCONTRIB_VERSION="1.049"
ARG JULIA_VERSION="1.9.3"
ARG YQ_VERSION="4.40.3"

# Set a label to keep track of all versions (these can be viewed with "docker inspect")
# Separate them with a '|' to potentially make them easily parsable

LABEL installed.software.versions="\
ROOT ${ROOT_VERSION}\
 | HepMC ${HEPMC_VERSION}\
 | Rivet ${RIVET_VERSION}\
 | Yoda ${YODA_VERSION}\
 | Fastjet ${FASTJET_VERSION}\
 | Fjcontrib ${FJCONTRIB_VERSION}\
 | Cppcheck ${CPPCHECK_VERSION}\
 | Cpplint ${CPPLINT_VERSION}\
 | Julia ${JULIA_VERSION}\
 | yq ${YQ_VERSION}"

# Add some info about authors and contents using the rules o the OCI Image Format Specification
# https://github.com/opencontainers/image-spec/blob/main/annotations.md#pre-defined-annotation-keys

LABEL org.opencontainers.image.title="SMASH-max"
LABEL org.opencontainers.image.authors="SMASH team"
LABEL org.opencontainers.image.url="https://smash-transport.github.io"
LABEL org.opencontainers.image.source="https://github.com/smash-transport/smash"
LABEL org.opencontainers.image.description="\
SMASH with ROOT ${ROOT_VERSION}, HepMC3 ${HEPMC_VERSION}, Rivet ${RIVET_VERSION} and several other tools"

# Update the distribution installing several additional packages non interactively

RUN <<EOF
  set -eux
  echo 'debconf debconf/frontend select Noninteractive' | debconf-set-selections
  apt-get update
  ln -sf /usr/share/zoneinfo/UTC /etc/localtime
  apt-get -y upgrade
  apt-get -y install \
    aspell  \
    aspell-en  \
    autoconf  \
    automake  \
    autotools-dev  \
    bc \
    build-essential  \
    bzip2  \
    ca-certificates  \
    clang-13 \
    clang-format-13 \
    clang-tidy-13 \
    clang-tools-13 \
    coreutils  \
    curl \
    cython3 \
    ddd  \
    debianutils  \
    diffutils  \
    doxygen  \
    evince  \
    findutils  \
    gawk  \
    g++ \
    gcc  \
    gcr  \
    gdb  \
    gfortran  \
    gnupg  \
    gnupg-utils  \
    gnuplot  \
    gnuplot-qt  \
    gperf  \
    graphviz \
    grep  \
    gzip  \
    imagemagick \
    ipython3  \
    jupyter-notebook \
    libc-bin  \
    libc-dev-bin  \
    libc-devtools \
    libc++-13-dev \
    libc++1-13 \
    libc++abi-13-dev \
    libc++abi1-13 \
    libclang-common-13-dev \
    libclang-cpp13 \
    libclang1-13 \
    libfreetype-dev \
    libgif-dev  \
    libgmp-dev  \
    libgraphviz-dev  \
    libjpeg-dev  \
    libjpeg8-dev  \
    libllvm13 \
    liblzma-dev  \
    libmpfr-dev  \
    libmpfrc++-dev  \
    libopenmpi-dev \
    libpcre3-dev \
    libpcre++-dev  \
    libpng-dev  \
    libpng-tools  \
    libssl-dev  \
    libtool  \
    libx11-dev  \
    libxext-dev \
    libxft-dev \
    libxml2-dev \
    libxpm-dev \
    libxxhash-dev \
    libzstd-dev \
    linux-libc-dev \
    llvm-13 \
    llvm-13-dev \
    llvm-13-linker-tools \
    llvm-13-runtime \
    llvm-13-tools \
    locales \
    locales-all \
    lsb-base  \
    lsb-release  \
    lzma  \
    python3-args  \
    python3-atomicwrites  \
    python3-colorama  \
    python3-cryptography  \
    python3-dateutil  \
    python3-decorator  \
    python3-dev  \
    python3-distro  \
    python3-distutils  \
    python3-entrypoints  \
    python3-future  \
    python3-importlib-metadata  \
    python3-ipython  \
    python3-ipython-genutils  \
    python-is-python3  \
    python3-mpi4py \
    python3-numexpr  \
    python3-numpy  \
    python3-pandas  \
    python3-patsy  \
    python3-pickleshare  \
    python3-pip \
    python3-pkg-resources  \
    python3-prompt-toolkit  \
    python3-psutil  \
    python3-ptyprocess  \
    python3-py  \
    python3-pytest  \
    python3-scipy  \
    python3-seaborn  \
    python3-secretstorage  \
    python3-setuptools  \
    python3-simplejson  \
    python3-statsmodels  \
    python3-systemd  \
    python3-tables  \
    python3-urllib3  \
    python3-yaml  \
    python3-zipp  \
    qiv  \
    sed  \
    sensible-utils  \
    sqlite3  \
    tar  \
    texlive-latex-base \
    texlive-pstricks \
    tree \
    unzip  \
    util-linux  \
    valgrind  \
    x11-common  \
    x11-xkb-utils  \
    x11-xserver-utils  \
    xauth  \
    xdot \
    xfonts-base  \
    xfonts-encodings  \
    xfonts-utils  \
    xkb-data  \
    xml-core  \
    xorg-sgml-doctools  \
    xserver-common  \
    xserver-xorg  \
    xserver-xorg-core  \
    xserver-xorg-input-all  \
    xserver-xorg-input-libinput  \
    xserver-xorg-legacy  \
    xxhash \
    xz-utils  \
    zenity  \
    zip
EOF

# Set PATH-like environment variables to comply with Docker suggestion to always define these
# before using them (otherwise an UndefinedVar warning is given by Docker at build time). Note
# that the base image already define LD_LIBRARY_PATH and we clearly do not want to reset this one.

ENV PYTHONPATH=''
ENV MANPATH=''
ENV CMAKE_PREFIX_PATH=''
ENV JUPYTER_PATH=''
ENV JUPYTER_CONFIG_DIR=''
ENV TEXMFHOME=''
ENV HOMETEXMF=''
ENV TEXMFCNF=''
ENV TEXINPUTS=''
ENV LATEXINPUTS=''

# Set cpu architecture flags for compilation
# ARG sets these flags as temporary environment variables that are used by cmake and by other tools

ARG CFLAGS="-O3 -march=${TARGET_ARCHITECTURE}"
ARG CPPFLAGS="${CFLAGS} -std=c++17"
ARG CXXFLAGS=${CPPFLAGS}

# Install ROOT

RUN <<EOF
  set -eux
  cd /SMASH
  wget https://root.cern/download/root_v${ROOT_VERSION}.source.tar.gz
  tar xf root_v${ROOT_VERSION}.source.tar.gz
  mv root-${ROOT_VERSION} root_src
  mkdir -p root_install root_build
  cd root_build
  cmake -DCMAKE_INSTALL_PREFIX=../root_install ../root_src -Dxrootd=OFF
  cmake --build . --target install -j$(nproc)
  printf '/SMASH/root_install/lib\n' >> /etc/ld.so.conf
  ldconfig
EOF

ENV ROOTSYS=/SMASH/root_install
ENV PATH=${ROOTSYS}/bin:${PATH}
ENV LD_LIBRARY_PATH=${ROOTSYS}/lib:${LD_LIBRARY_PATH}
ENV SHLIB_PATH=${ROOTSYS}/lib
ENV LIBPATH=${ROOTSYS}/lib
ENV PYTHONPATH=${ROOTSYS}/lib:${PYTHONPATH}
ENV MANPATH=${ROOTSYS}/man:${MANPATH}
ENV CMAKE_PREFIX_PATH=${ROOTSYS}:${CMAKE_PREFIX_PATH}
ENV JUPYTER_PATH=${ROOTSYS}/etc/notebook:${JUPYTER_PATH}
ENV JUPYTER_CONFIG_DIR=${ROOTSYS}/etc/notebook:${JUPYTER_CONFIG_DIR}
ENV CLING_STANDARD_PCH=none

# Install HepMC3

ARG HEPMC_INSTALL="/SMASH/hepmc3_install"

RUN <<EOF
  set -eux
  cd /SMASH
  wget https://hepmc.web.cern.ch/hepmc/releases/HepMC3-${HEPMC_VERSION}.tar.gz
  tar xf HepMC3-${HEPMC_VERSION}.tar.gz
  mkdir -p HepMC3_build
  cd HepMC3_build
  cmake -DHEPMC3_ENABLE_ROOTIO:BOOL=ON \
        -DHEPMC3_ENABLE_TEST:BOOL=OFF \
        -DHEPMC3_INSTALL_INTERFACES:BOOL=ON \
        -DHEPMC3_ENABLE_PYTHON:BOOL=OFF \
        -DHEPMC3_BUILD_STATIC_LIBS:BOOL=OFF \
        -DHEPMC3_BUILD_DOCS:BOOL=OFF \
        -DCMAKE_INSTALL_PREFIX=../hepmc3_install ../HepMC3-${HEPMC_VERSION}
  make -j$(nproc)
  make install
  printf "${HEPMC_INSTALL}/lib\n" >> /etc/ld.so.conf
  ldconfig
EOF

ENV PATH=${HEPMC_INSTALL}/bin:${PATH}
# This is needed to make ROOT correctly find HEPMC files
ENV ROOT_INCLUDE_PATH=/SMASH/hepmc3_install/include/

# Install gitstats

RUN <<EOF
  set -eux
  cd /SMASH
  git clone https://github.com/hoxu/gitstats.git
  cd gitstats
  make install
EOF

# Install cppcheck and cpplint

RUN <<EOF
  set -eux
  cd /SMASH
  wget https://github.com/danmar/cppcheck/archive/${CPPCHECK_VERSION}.tar.gz
  tar xf ${CPPCHECK_VERSION}.tar.gz
  cd cppcheck-${CPPCHECK_VERSION}
  mkdir -p build
  cd build
  cmake ..
  cmake --build .
  cd /SMASH
  wget https://github.com/cpplint/cpplint/archive/refs/tags/${CPPLINT_VERSION}.tar.gz
  tar xf ${CPPLINT_VERSION}.tar.gz
  mkdir -p /SMASH/auxiliary_progs/
  mv cppcheck-${CPPCHECK_VERSION}/build/bin/cppcheck cpplint-${CPPLINT_VERSION}/cpplint.py /SMASH/auxiliary_progs/
  chmod +x /SMASH/auxiliary_progs/*
EOF

ENV PATH=${PATH}:/SMASH/auxiliary_progs

# Install Rivet (the versions of the dependent libraries listed at the top of this file
# are chosen according to the bootstrap script downloadable from the Rivet website)

ARG RIVET_INSTALL_PREFIX="/SMASH/Rivet_install"
ARG RIVET_BUILD_PREFIX="/SMASH/rivet_build"
ARG FASTJET_PATH=${RIVET_INSTALL_PREFIX}
ARG YODA_PATH=${RIVET_INSTALL_PREFIX}

RUN <<EOF
  set -eux -o pipefail
  printf 'Preparing Rivet installation\n'
  mkdir -p ${RIVET_INSTALL_PREFIX}/etc/bash_completion.d \
           ${RIVET_INSTALL_PREFIX}/lib
  ln -s lib ${RIVET_INSTALL_PREFIX}/lib64
  mkdir -p ${RIVET_BUILD_PREFIX}
  printf 'Installing Fastjet\n'
  cd ${RIVET_BUILD_PREFIX}
  wget --no-check-certificate --no-verbose \
       http://www.fastjet.fr/repo/fastjet-${FASTJET_VERSION}.tar.gz -O- | tar --no-same-owner -xz
  cd fastjet-${FASTJET_VERSION}
  ./configure --prefix=${RIVET_INSTALL_PREFIX} --enable-shared --disable-auto-ptr --enable-allcxxplugins
  make -j$(nproc)
  make install
  printf 'Installing Fastjet-contrib\n'
  cd ${RIVET_BUILD_PREFIX}
  wget --no-check-certificate --no-verbose \
       http://cedar-tools.web.cern.ch/downloads/fjcontrib-${FJCONTRIB_VERSION}.tar.gz -O- | tar --no-same-owner -xz
  cd fjcontrib-${FJCONTRIB_VERSION}
  ./configure --fastjet-config=${FASTJET_PATH}/bin/fastjet-config CXXFLAGS=-fPIC
  sed -i.bak -e 's/#include <array>/#include <array>\n#include <limits>/' LundPlane/EEHelpers.hh
  make -j1 fragile-shared-install
  make install
  printf 'Installing Yoda\n'
  cd ${RIVET_BUILD_PREFIX}
  wget --no-check-certificate --no-verbose \
       http://cedar-tools.web.cern.ch/downloads/YODA-${YODA_VERSION}.tar.gz -O- | tar --no-same-owner -xz
  cd YODA-${YODA_VERSION}
  ./configure --prefix=${YODA_PATH}
  make -j$(nproc)
  make install
  cp yodaenv.sh ${YODA_PATH}/yodaenv.sh
  printf 'Installing Rivet\n'
  cd ${RIVET_BUILD_PREFIX}
  wget --no-check-certificate --no-verbose \
       http://cedar-tools.web.cern.ch/downloads/Rivet-${RIVET_VERSION}.tar.gz -O- | tar --no-same-owner -xz
  cd Rivet-${RIVET_VERSION}
  ./configure --prefix=${RIVET_INSTALL_PREFIX} \
              --with-yoda=${YODA_PATH} \
              --with-hepmc3=${HEPMC_INSTALL} \
              --with-fastjet=${FASTJET_PATH}
  make -j$(nproc)
  make install
  cp rivetenv.sh ${RIVET_INSTALL_PREFIX}
  printf "${RIVET_INSTALL_PREFIX}/lib\n" >> /etc/ld.so.conf
  ldconfig
EOF

ENV PYTHON=/usr/bin/python3
ENV PATH=${RIVET_INSTALL_PREFIX}/bin:${PATH}
ENV LD_LIBRARY_PATH="${RIVET_INSTALL_PREFIX}/lib:${RIVET_INSTALL_PREFIX}/lib64:${LD_LIBRARY_PATH}"
# The python version must be updated together with the underlying linux distribution.
# Unfortunately something like $(python3 -V | cut -d " " -f 2) does not work in an ENV statement.
ENV PYTHONPATH="${RIVET_INSTALL_PREFIX}/lib/python3.8.10/site-packages:${PYTHONPATH}"
ENV TEXMFHOME="${RIVET_INSTALL_PREFIX}/share/Rivet/texmf:${TEXMFHOME}"
ENV HOMETEXMF="${RIVET_INSTALL_PREFIX}/share/Rivet/texmf:${HOMETEXMF}"
ENV TEXMFCNF="${RIVET_INSTALL_PREFIX}/share/Rivet/texmf/cnf:${TEXMFCNF}"
ENV TEXINPUTS="${RIVET_INSTALL_PREFIX}/share/Rivet/texmf/tex:${TEXINPUTS}"
ENV LATEXINPUTS="${RIVET_INSTALL_PREFIX}/share/Rivet/texmf/tex:${LATEXINPUTS}"

# Install julia

ARG JULIA_PATH=/opt/julia-${JULIA_VERSION}

RUN <<EOF
  set -eux
  cd /opt
  wget https://github.com/JuliaLang/julia/releases/download/v${JULIA_VERSION}/julia-${JULIA_VERSION}.tar.gz
  tar xf julia-${JULIA_VERSION}.tar.gz
  cd julia-${JULIA_VERSION}
  make -j$(nproc)
EOF

ENV PATH=${JULIA_PATH}:${PATH}

# Install yq

ARG YQ_BINARY=yq_linux_amd64

RUN <<EOF
  set -eux
  cd /opt
  wget https://github.com/mikefarah/yq/releases/download/v${YQ_VERSION}/${YQ_BINARY} -O /usr/bin/yq
  chmod +x /usr/bin/yq
EOF

# Recompile SMASH (if the repository exists)

RUN <<EOF
  set -eux
  if [[ -d /SMASH/smash ]]
  then
    cd /SMASH
    rm -rf smash_bin
    mkdir -p smash_bin
    cd smash_bin
    cmake -DPythia_CONFIG_EXECUTABLE=/SMASH/pythia_install/bin/pythia8-config\
          -DCMAKE_INSTALL_PREFIX=${HEPMC_INSTALL} /SMASH/smash
    make -j$(nproc)
    cd /SMASH
    if [ ! -f /SMASH/smash_bin/smash ]
    then
      rm -rf /SMASH/smash /SMASH/smash_bin
    fi
  fi
EOF

# Set the locations of the SMASH sources and build directory for the variables used by FindSMASH.cmake

ENV SMASH_DIR=/SMASH/smash
ENV SMASH_BUILD_DIR=/SMASH/smash_bin

# Miscellaneous:
#  1. Remove most of the sources
#  2. Allow Imagemagick to process PS, EPS and PDF files in case the user needs to use `convert` inside the container
#  3. Set a symbolic link to clang-format-13 and clang-tidy-13 to make them usable from the SMASH formatting script
#  4. Install cmake-format

RUN <<EOF
  set -eux
  cd /SMASH
  rm -rf HepMC3-${HEPMC_VERSION} \
         HepMC3-${HEPMC_VERSION}.tar.gz \
         HepMC3_build \
         gitstats \
         cppcheck-${CPPCHECK_VERSION} \
         ${CPPCHECK_VERSION}.tar.gz \
         cpplint-${CPPLINT_VERSION} \
         ${CPPLINT_VERSION}.tar.gz \
         root_v${ROOT_VERSION}.source.tar.gz \
         root_build \
         root_src \
         ${RIVET_BUILD_PREFIX} \
         /var/cache/apt/archives/*.deb
  sed -i -e "/ghostscript/d" -e "/PS/d" -e "/PS2/d" -e "/PS3/d" \
         -e "/EPS/d" -e "/PDF/d" -e "/XPS/d" /etc/ImageMagick-6/policy.xml
  ln -s /usr/bin/clang-format-13 /usr/bin/clang-format
  ln -s /usr/bin/clang-tidy-13 /usr/bin/clang-tidy
  pip install cmakelang
EOF

# Set a few other environment variables (just to be sure which locale is used)

ENV LC_ALL=en_US.UTF-8
ENV LANG en=US.UTF-8
ENV LANGUAGE=en_US.UTF-8
ENV PROMPT_CONTAINER_NAME='SMASH-max'
