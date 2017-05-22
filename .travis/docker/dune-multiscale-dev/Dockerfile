FROM dunecommunity/dailywork_debian-full-clion

MAINTAINER René Milk <rene.milk@wwu.de>

#openssh server needed for orte daemon
RUN apt-get update && \
    apt-get remove -y libeigen3-dev libscotchmetis-dev libscotchparmetis-dev && \
    apt-get autoremove -y && \
    apt-get install -yq --no-install-recommends cmake-curses-gui \
            libfftw3-mpi-dev libfftw3-bin libfftw3-dev libopenmpi-dev \
            libglpk-dev libfreetype6 mpi-default-dev \
            python-pip python3-pip python3-venv python-virtualenv \
            openssh-server

ENV PATH /usr/lib/ccache:${PATH} \
    SRC_DCTRL ${SUPERDIR}/dune-common/bin/dunecontrol
