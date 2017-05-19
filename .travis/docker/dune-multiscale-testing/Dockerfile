FROM dunecommunity/testing-base

# COMPILER Selection (ARGS cannot be uppercase)
ARG cc
ARG cxx
ENV CC $cc
ENV CXX $cxx

ENV MY_MODULE dune-multiscale
ENV OPTS config.opts/travis.ninja
ENV SRC_DCTRL ${SUPERDIR}/dune-common/bin/dunecontrol

WORKDIR ${SUPERDIR}
ADD . ${SUPERDIR}


RUN apt-get update && \
    apt-get remove -y libeigen3-dev libscotchmetis-dev libscotchparmetis-dev && \
    apt-get autoremove -y && \
    apt-get install -yq --no-install-recommends cmake-curses-gui \
            libfftw3-mpi-dev libfftw3-bin libfftw3-dev libopenmpi-dev \
            libglpk-dev libfreetype6 mpi-default-dev \
            python-pip python3-pip python3-venv python-virtualenv && \
    cd ${SUPERDIR} && \
    git submodule foreach git reset --hard && \
    git submodule update --init --recursive && \
    ./local/bin/download_external_libraries.py && \
    ./local/bin/build_external_libraries.py  && \
    ${SRC_DCTRL} ${BLD} all && \
    rm -rf ${MY_MODULE} && \
    mkdir ${MY_MODULE} && \
    rm -rf ${DUNE_BUILD_DIR}/${MY_MODULE}/CMakeCache.txt
# THE END
ENV DEBIAN_FRONTEND teletype
