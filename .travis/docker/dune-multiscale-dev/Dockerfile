FROM dunecommunity/dailywork_debian-full-clion

MAINTAINER René Milk <rene.milk@wwu.de>

RUN apt-get update && \
    apt-get remove -y libeigen3-dev libscotchmetis-dev libscotchparmetis-dev && \
    apt-get autoremove -y && \
    apt-get install -yq --no-install-recommends cmake-curses-gui \
            libfftw3-mpi-dev libfftw3-bin libfftw3-dev libopenmpi-dev \
            libglpk-dev libfreetype6 mpi-default-dev \
            python-pip python3-pip python3-venv python-virtualenv

# COMPILER Selection (ARGS cannot be uppercase)
ARG cc
ARG cxx
ENV CC $cc \
    CXX $cxx

ENV MY_MODULE dune-multiscale \
    PATH /usr/lib/ccache:${PATH} \
    OPTS config.opts/stinson/${CC}_debug \
    SRC_DCTRL ${SUPERDIR}/dune-common/bin/dunecontrol \
    DSETUP dune_setup ${OPTS}
