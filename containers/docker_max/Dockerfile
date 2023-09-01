# syntax=docker/dockerfile:1.4

# Import the base smash Docker container image.
# Adjust the address as needed. In particular, be sure to use a base image compatible with your cpu.
# The official Docker images are prepared for a generic Intel/AMD x86-64 cpu.

FROM ghcr.io/smash-transport/smash:newest

# Set bash as the default shell in this Dockerfile.

SHELL ["/bin/bash", "-c"]

# Set the proper general cpu family architecture flag.
# Please, for more information refer to the documention of your compiler.

ARG TARGET_ARCHITECTURE="x86-64"

# Add some info about authors and contents.

LABEL maintainer.name="SMASH team"
LABEL maintainer.website="https://smash-transport.github.io"
LABEL description="SMASH with ROOT 6.26.10, HepMC3 3.2.6, Rivet, 3.1.7 and several other tools"
LABEL org.opencontainers.image.source https://github.com/smash-transport/smash

# Install the packages non interactively.

RUN echo 'debconf debconf/frontend select Noninteractive' | debconf-set-selections

# Update the distribution and install several additional packages.

RUN apt-get update && \
ln -sf /usr/share/zoneinfo/UTC /etc/localtime && \
apt-get -y upgrade && \
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
screen  \
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

# Set cpu architecture flags for compilation
# ARG sets these flags as temporary environment variables that are used by cmake and by other tools

ARG CFLAGS="-O3 -march=${TARGET_ARCHITECTURE}"
ARG CPPFLAGS="${CFLAGS} -std=c++17"
ARG CXXFLAGS=${CPPFLAGS}

# Install ROOT.

ARG rootV="6.26.10"

RUN cd /SMASH && \
wget https://root.cern/download/root_v${rootV}.source.tar.gz && \
tar xf root_v${rootV}.source.tar.gz && \
mv root-${rootV} root_src && \
mkdir -p root_install && \
mkdir -p root_build && \
cd root_build && \
cmake -DCMAKE_INSTALL_PREFIX=../root_install ../root_src -Dxrootd=OFF && \
cmake --build . --target install -j$(nproc)
ENV ROOTSYS=/SMASH/root_install
ENV PATH=$ROOTSYS/bin:$PATH
ENV LD_LIBRARY_PATH=$ROOTSYS/lib:$LD_LIBRARY_PATH
ENV SHLIB_PATH=$ROOTSYS/lib
ENV LIBPATH=$ROOTSYS/lib
ENV PYTHONPATH=$ROOTSYS/lib:$PYTHONPATH
ENV MANPATH=$ROOTSYS/man:$MANPATH
ENV CMAKE_PREFIX_PATH=$ROOTSYS:$CMAKE_PREFIX_PATH
ENV JUPYTER_PATH=$ROOTSYS/etc/notebook:$JUPYTER_PATH
ENV JUPYTER_CONFIG_DIR=$ROOTSYS/etc/notebook:$JUPYTER_CONFIG_DIR
ENV CLING_STANDARD_PCH=none
RUN echo /SMASH/root_install/lib >> /etc/ld.so.conf && ldconfig

# Install HepMC3.

ARG hepmcV="3.2.6"

RUN cd /SMASH && \
wget https://hepmc.web.cern.ch/hepmc/releases/HepMC3-${hepmcV}.tar.gz && \
tar xf HepMC3-${hepmcV}.tar.gz && \
mkdir -p HepMC3_build && \
cd HepMC3_build && \
cmake -DHEPMC3_ENABLE_ROOTIO:BOOL=ON -DHEPMC3_ENABLE_TEST:BOOL=OFF -DHEPMC3_INSTALL_INTERFACES:BOOL=ON -DHEPMC3_ENABLE_PYTHON:BOOL=OFF -DHEPMC3_BUILD_STATIC_LIBS:BOOL=OFF -DHEPMC3_BUILD_DOCS:BOOL=OFF -DCMAKE_INSTALL_PREFIX=../hepmc3_install ../HepMC3-${hepmcV} && \
make -j$(nproc) && \
make install
ENV PATH="/SMASH/hepmc3_install/bin:$PATH"
RUN echo /SMASH/hepmc3_install/lib >> /etc/ld.so.conf && ldconfig

# Install gitstats.

RUN cd /SMASH && \
git clone https://github.com/hoxu/gitstats.git && \
cd gitstats && \
make install

# Install cppcheck and cpplint.

ARG cppcheckV="2.8"
ARG cpplintV="1.6.0"

RUN cd /SMASH && \
wget https://github.com/danmar/cppcheck/archive/${cppcheckV}.tar.gz && \
tar xf ${cppcheckV}.tar.gz && \
cd cppcheck-${cppcheckV} && \
mkdir -p build && cd build && cmake .. && cmake --build . && \
cd /SMASH && wget https://github.com/cpplint/cpplint/archive/refs/tags/${cpplintV}.tar.gz && tar xf ${cpplintV}.tar.gz && \
mkdir -p /SMASH/auxiliary_progs/ && mv cppcheck-${cppcheckV}/build/bin/cppcheck cpplint-${cpplintV}/cpplint.py /SMASH/auxiliary_progs/  && \
chmod +x /SMASH/auxiliary_progs/*
ENV PATH=$PATH:/SMASH/auxiliary_progs

# Install Rivet.
# The versions of the dependent libraries are chosen according to the bootstrap script downloadable from the Rivet website.

ARG RIVET_VERSION="3.1.7"
ARG YODA_VERSION="1.9.7"
ARG FASTJET_VERSION="3.4.0"
ARG FJCONTRIB_VERSION="1.049"
ARG CXXDIALECT="c++17"
ARG INSTALL_PREFIX="/SMASH/Rivet_install"
ARG HEPMCPATH="/SMASH/hepmc3_install"
ARG BUILD_PREFIX="/SMASH/rivet_build"

ENV PYTHON=/usr/bin/python3
ENV PATH=/SMASH/Rivet_install/bin:$PATH
ENV LD_LIBRARY_PATH="/SMASH/Rivet_install/lib:/SMASH/Rivet_install/lib64:$LD_LIBRARY_PATH"

RUN mkdir -p $INSTALL_PREFIX/etc/bash_completion.d && mkdir -p $INSTALL_PREFIX/lib && ln -s lib $INSTALL_PREFIX/lib64

RUN mkdir -p $BUILD_PREFIX && cd $BUILD_PREFIX

RUN wget --no-check-certificate --no-verbose http://www.fastjet.fr/repo/fastjet-$FASTJET_VERSION.tar.gz -O- | tar --no-same-owner -xz &&\
cd fastjet-$FASTJET_VERSION && ./configure --prefix=$INSTALL_PREFIX --enable-shared --disable-auto-ptr --enable-allcxxplugins &&\
make -j$(nproc) && make install
ARG FASTJETPATH=$INSTALL_PREFIX

RUN cd $BUILD_PREFIX &&\
wget --no-check-certificate --no-verbose http://cedar-tools.web.cern.ch/downloads/fjcontrib-$FJCONTRIB_VERSION.tar.gz -O- | tar --no-same-owner -xz &&\
cd fjcontrib-$FJCONTRIB_VERSION && ./configure --fastjet-config=$FASTJETPATH/bin/fastjet-config CXXFLAGS=-fPIC &&\
sed -i.bak -e 's/#include <array>/#include <array>\n#include <limits>/' LundPlane/EEHelpers.hh &&\
make -j1 fragile-shared-install && make install

RUN cd $BUILD_PREFIX &&\
wget --no-check-certificate --no-verbose http://cedar-tools.web.cern.ch/downloads/YODA-$YODA_VERSION.tar.gz -O- | tar --no-same-owner -xz &&\
cd YODA-$YODA_VERSION && ./configure --prefix=$INSTALL_PREFIX && make -j$(nproc) && make install && cp yodaenv.sh $INSTALL_PREFIX/yodaenv.sh
ARG YODAPATH=$INSTALL_PREFIX

# Remove all the experimental analyses and leave only the analysis for simulators.
# The lines need to be checked and probably adjusted at each new Rivet release.

RUN cd $BUILD_PREFIX &&\
wget --no-check-certificate --no-verbose http://cedar-tools.web.cern.ch/downloads/Rivet-$RIVET_VERSION.tar.gz -O- | tar --no-same-owner -xz &&\
cd Rivet-$RIVET_VERSION &&\
./configure --prefix=$INSTALL_PREFIX --with-yoda=$YODAPATH -with-hepmc3=$HEPMCPATH --with-fastjet=$FASTJETPATH &&\
make -j$(nproc) && make install && cp rivetenv.sh $INSTALL_PREFIX

RUN echo $INSTALL_PREFIX/lib >> /etc/ld.so.conf && ldconfig 

# The python version must be updated together with the underlying linux distribution.
# Unfortunately something like $(python3 -V | cut -d " " -f 2) does not work in an ENV statement.

ENV PYTHONPATH="$INSTALL_PREFIX/lib/python3.8.10/site-packages:$PYTHONPATH"
ENV TEXMFHOME="$INSTALL_PREFIX/share/Rivet/texmf:$TEXMFHOME"
ENV HOMETEXMF="$INSTALL_PREFIX/share/Rivet/texmf:$HOMETEXMF"
ENV TEXMFCNF="$INSTALL_PREFIX/share/Rivet/texmf/cnf:$TEXMFCNF"
ENV TEXINPUTS="$INSTALL_PREFIX/share/Rivet/texmf/tex:$TEXINPUTS"
ENV LATEXINPUTS="$INSTALL_PREFIX/share/Rivet/texmf/tex:$LATEXINPUTS"

# Install julia.

ARG julia_majV="1.8"
ARG julia_minV="5"

RUN cd /opt &&\
wget https://github.com/JuliaLang/julia/releases/download/v${julia_majV}.${julia_minV}/julia-${julia_majV}.${julia_minV}.tar.gz &&\
tar xf julia-${julia_majV}.${julia_minV}.tar.gz && cd julia-${julia_majV}.${julia_minV} && make -j$(nproc)
ARG julia_path=/opt/julia-${julia_majV}.${julia_minV}
ENV PATH=${julia_path}:$PATH

# Recompile SMASH (if the repository exists).

RUN <<EOF
if [ -d /SMASH/smash ]
then
  cd /SMASH
  rm -rf smash_bin
  mkdir -p smash_bin
  cd smash_bin
  cmake /SMASH/smash -DPythia_CONFIG_EXECUTABLE=/SMASH/pythia_install/bin/pythia8-config -DCMAKE_INSTALL_PREFIX=/SMASH/hepmc3_install \
  || echo "cmake SMASH build configuration failed!"
  make -j$(nproc) || echo "SMASH or tests compilation failed!"
  cd /SMASH
  if [ ! -f /SMASH/smash_bin/smash ]
  then
    rm -rf /SMASH/smash /SMASH/smash_bin
  fi
fi
EOF

# Set the locations of the SMASH sources and build directory for the variables used by FindSMASH.cmake.

ENV SMASH_DIR=/SMASH/smash
ENV SMASH_BUILD_DIR=/SMASH/smash_bin

# Remove most of the sources.

RUN cd /SMASH && \
rm -rf HepMC3-${hepmcV} HepMC3-${hepmcV}.tar.gz HepMC3_build pythia${pythiaV} gitstats cppcheck-${cppcheckV} ${cppcheckV}.tar.gz \
cpplint-${cpplintV} ${cpplintV}.tar.gz root_v${rootV}.source.tar.gz root_build root_src $BUILD_PREFIX /var/cache/apt/archives/*.deb

# Allow Imagemagick to process PS, EPS and PDF files in case the user needs to use `convert` inside the container

RUN sed -i -e "/ghostscript/d" -e "/PS/d" -e "/PS2/d" -e "/PS3/d" -e "/EPS/d" -e "/PDF/d" -e "/XPS/d" /etc/ImageMagick-6/policy.xml

# Set a few other environment variables (just to be sure which locale is used).

ENV LC_ALL=en_US.UTF-8
ENV LANG en=US.UTF-8
ENV LANGUAGE=en_US.UTF-8

# We set a symbolic link to clang-format-13 and clang-tidy-13 to make them usable from the SMASH formatting script

RUN ln -s /usr/bin/clang-format-13 /usr/bin/clang-format &&\
ln -s /usr/bin/clang-tidy-13 /usr/bin/clang-tidy

# we install cmake-format

RUN pip install cmakelang

# Shell prompt when using the container image under Docker.

ENV PS1="S - \[\e]0;\u@\h: \w\a\]${debian_chroot:+($debian_chroot)}\u@\h:\w\$ "

# Shell prompt when using the container image under Singularity (CUSTOM_ENV file automatically sourced by Singularity).

RUN mkdir -p /.singularity.d/env/ && \
CUSTOM_ENV=/.singularity.d/env/99-zz_custom_env.sh && \
echo '#!/bin/bash' > $CUSTOM_ENV && \
echo 'PS1="S - \[\e]0;\u@\h: \w\a\]${debian_chroot:+($debian_chroot)}\u@\h:\w\$ "' >> $CUSTOM_ENV && \
chmod 755 $CUSTOM_ENV

# Build wrapper script around screen executable to make it usable from within the container
# (it needs write access to SCREENDIR folder)

RUN mv /usr/bin/screen /usr/bin/screen.exe && \
cat <<EOF > /usr/bin/screen
printf '%s\n'\\
       'Executing screen wrapper script around original program /usr/bin/screen.exe'\\
       'We try to set the environment variable SCREENDIR to \${HOME}/.screen_dir,'\\
       'if this directory is not accessible, please, manually set SCREENDIR to another directory'\\
       'before launching /usr/bin/screen.exe'
sleep 3
export SCREENDIR="\${HOME}/.screen_dir"
mkdir -p "\${SCREENDIR}"
chmod 700 "\${SCREENDIR}"
/usr/bin/screen.exe "\$@"
EOF
RUN chmod 755 /usr/bin/screen
