# syntax=docker/dockerfile:1.4

FROM ubuntu:22.04

# Set bash as the default shell in this Dockerfile

SHELL ["/bin/bash", "-c"]

# Set the proper general cpu family architecture flag.
# Please, for more information refer to the documentation of your compiler.
# The default value is set to "x86-64-v2" which is a conservative choice for modern
# x86-64 CPUs (cf. https://en.wikipedia.org/wiki/X86-64#Microarchitecture_levels)

ARG TARGET_ARCHITECTURE="x86-64-v2"

# Add some info about authors and contents using the rules o the OCI Image Format Specification
# https://github.com/opencontainers/image-spec/blob/main/annotations.md#pre-defined-annotation-keys

LABEL org.opencontainers.image.title="SMASH-mini"
LABEL org.opencontainers.image.authors="SMASH team"
LABEL org.opencontainers.image.url="https://smash-transport.github.io"
LABEL org.opencontainers.image.source="https://github.com/smash-transport/smash"
LABEL org.opencontainers.image.description="\
Basic SMASH image tuned for ${TARGET_ARCHITECTURE} architecture\
 (sources, compiler and bare executable, without ROOT, HepMC3 and Rivet)."

# Fail if TARGET_ARCHITECTURE is set to "native" (this is not allowed in this Dockerfile
# since it would make the image non portable)

RUN <<EOF
  if [[ "$TARGET_ARCHITECTURE" = "native" ]]; then
    echo "Warning: TARGET_ARCHITECTURE is set to 'native' - be sure this is desired." >&2
  fi
EOF

# Install minimal set of requirements

RUN <<EOF
  set -eux
  apt-get update
  apt-get -y upgrade
  DEBIAN_FRONTEND=noninteractive apt-get -y install \
    bash \
    cmake \
    emacs \
    git \
    g++ \
    libeigen3-dev \
    libgsl-dev \
    rsync \
    screen \
    wget \
    vim
EOF

ARG C_CXX_ARCH_FLAG="-march=${TARGET_ARCHITECTURE}"

# Set LD_LIBRARY_PATH environment variable to comply with Docker suggestion to always define
# variables before using them (otherwise an UndefinedVar warning is given by Docker at build time).
# Note that the base Ubuntu image already define PATH and we clearly do not want to reset this one.

ENV LD_LIBRARY_PATH=''

# Create main project directory

WORKDIR /SMASH

# Create directory for mounting filesystems

RUN mkdir /MNT

# Get and build public SMASH

ARG SMASH_TAG="latest"
ENV PYTHIA8=/SMASH/pythia_install/

# Select between cloning SMASH publicly or using local development sources

ARG USE_LOCAL_SOURCE=false

# syntax=docker/dockerfile:1

ARG SMASH_TAG=latest
ARG USE_LOCAL_SOURCE=false

ENV PYTHIA8=/SMASH/pythia_install/
# syntax=docker/dockerfile:1

ENV PYTHIA8=/SMASH/pythia_install/

RUN --mount=from=local_smash,target=/__local_smash <<'EOF'
  set -euxo pipefail
  mkdir -p /SMASH
  cd /SMASH

  if [[ "${USE_LOCAL_SOURCE}" == "true" ]]; then
    echo "Using local SMASH development sources..."
    cp -a /__local_smash /SMASH/smash
  else
    echo "Cloning public SMASH repository..."
    git clone https://github.com/smash-transport/smash.git smash
    cd smash

    if [[ "${SMASH_TAG}" == "latest" ]]; then
      latest_tag=$(git describe --tags "$(git rev-list --tags --max-count=1)")
      git switch --detach "${latest_tag}"
    else
      git checkout "${SMASH_TAG}"
    fi
  fi

  cd /SMASH/smash

  find_pythia_strings=($(grep 'find_package(Pythia' src/CMakeLists.txt | grep -Eo '[1-9]?[0-9]\.[0-9]{3}'))

  if [[ ${#find_pythia_strings[@]} -ne 1 ]]; then
    echo "Error: Exactly one Pythia version expected from src/CMakeLists.txt, found: ${#find_pythia_strings[@]}" >&2
    exit 1
  fi

  PYTHIA_VERSION="${find_pythia_strings[0]}"
  PYTHIA_VERSION_NODOT="${PYTHIA_VERSION/./}"

  if [[ -z "${PYTHIA_VERSION}" ]]; then
    echo "Error: could not detect Pythia version from src/CMakeLists.txt" >&2
    exit 1
  fi

  echo "Detected Pythia version: ${PYTHIA_VERSION}"
  echo "Using tarball version: ${PYTHIA_VERSION_NODOT}"

  cd /SMASH
  wget "https://pythia.org/download/pythia83/pythia${PYTHIA_VERSION_NODOT}.tgz"
  tar xf "pythia${PYTHIA_VERSION_NODOT}.tgz"
  rm "pythia${PYTHIA_VERSION_NODOT}.tgz"

  cd "/SMASH/pythia${PYTHIA_VERSION_NODOT}"
  ./configure --cxx-common="${C_CXX_ARCH_FLAG} -O3 -std=c++17 -fPIC -pthread" --prefix="${PYTHIA8}"
  make -j"$(nproc)"
  make install

  echo "${PYTHIA8}/lib" > /etc/ld.so.conf.d/pythia8.conf
  ldconfig

  mkdir -p /SMASH/smash/build
  mkdir -p /SMASH/smash_bin
  cd /SMASH/smash/build

  cmake \
    -DPythia_CONFIG_EXECUTABLE="/SMASH/pythia${PYTHIA_VERSION_NODOT}/bin/pythia8-config" \
    -DSMASH_ARCH_FLAG="${C_CXX_ARCH_FLAG}" \
    ..

  make smash -j"$(nproc)"

  if [[ -f /SMASH/smash/build/smash ]]; then
    cp /SMASH/smash/build/smash /SMASH/smash_bin/
    cp /SMASH/smash/build/config.yaml /SMASH/smash_bin/
    rm -rf /SMASH/smash/build
  else
    echo "SMASH executable not found! Something went wrong compiling SMASH!" >&2
    exit 1
  fi
EOF

ENV PYTHIA8DATA=${PYTHIA8}/share/Pythia8/xmldoc
ENV PYTHIA_LIBRARY=${PYTHIA8}/lib
ENV PYTHIA_INCLUDE=${PYTHIA8}/include
ENV LD_LIBRARY_PATH=${PYTHIA_LIBRARY}:${LD_LIBRARY_PATH}
ENV PATH=${PYTHIA8}/bin:/SMASH/smash_bin:${PATH}



# Singularity (CUSTOM_ENV file automatically sourced by Singularity).

RUN <<'EOF'
  set -eux
  COMMON_PS1='${PROMPT_CONTAINER_NAME:-SMASH} - \D{%d.%m.%y %H:%M:%S} \001\033[1;94m\002\w\001\033[0m\002 \$ '
  printf "\n\nexport PS1='🐳 %s'\n" "${COMMON_PS1}" >> ~/.bashrc
  mkdir -p /.singularity.d/env
  CUSTOM_ENV=/.singularity.d/env/99-zz_custom_env.sh
  printf "%s\n\n\nexport PS1='📦 %s'\n\n"\
         '#!/bin/bash'\
         "${COMMON_PS1}" > "${CUSTOM_ENV}"
  chmod 755 "${CUSTOM_ENV}"
EOF

ENV PROMPT_CONTAINER_NAME='SMASH-min'

# Build wrapper script around screen executable to make it usable from within the container
# (it needs write access to SCREENDIR folder). Note that quoting the sentinel word in the
# here-doc construct prevents to expand any ${...} inside (which is what we want here).

RUN <<EOF
  set -eux
  mv /usr/bin/screen /usr/bin/screen.exe
  cat > /usr/bin/screen <<'INTERNAL_EOF'
    printf '%s\n'\
          'Executing screen wrapper script around original program /usr/bin/screen.exe'\
          'We try to set the environment variable SCREENDIR to ${HOME}/.screen_dir,'\
          'if this directory is not accessible, please, manually set SCREENDIR to another directory'\
          'before launching /usr/bin/screen.exe'
    sleep 3
    export SCREENDIR="${HOME}/.screen_dir"
    mkdir -p "${SCREENDIR}"
    chmod 700 "${SCREENDIR}"
    /usr/bin/screen.exe "$@"
INTERNAL_EOF
  chmod 755 /usr/bin/screen
EOF
