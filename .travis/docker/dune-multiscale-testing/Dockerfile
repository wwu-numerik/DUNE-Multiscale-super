ARG BASE
FROM dunecommunity/testing-base_${BASE}

ENV MY_MODULE dune-multiscale

ARG COMMIT

RUN git clone https://github.com/wwu-numerik/dune-multiscale-super.git ${SUPERDIR} && \
    cd ${SUPERDIR} && \
    git checkout ${COMMIT} && \
    git submodule foreach git reset --hard && \
    git submodule update --init --recursive && \
    git submodule status
WORKDIR ${SUPERDIR}

RUN apt update && apt install -y libfftw3-mpi-dev libfftw3-mpi3 libfftw3-dev || echo "on arch"
RUN ${PACMAN_UPDATE} && fftw || echo "on debian"

# TRAVIS before_install:
ENV OPTS config.opts/travis.ninja

ENV PATH /usr/lib/ccache:${PATH}
# our local scripts look for an OPTS env entry
RUN ./bin/download_external_libraries.py
RUN ./bin/build_external_libraries.py
RUN virtualenv --python=python3 ./setup-venv && . ./setup-venv/bin/activate

# TRAVIS install
ENV SRC_DCTRL ${SUPERDIR}/dune-common/bin/dunecontrol
RUN ${SRC_DCTRL} ${BLD} all

# ensures ${MY_MODULE} from travis own checkout is used via clean mount point, but enables better compile re-use
RUN rm -rf ${MY_MODULE}
RUN mkdir ${MY_MODULE}
RUN rm -rf ${DUNE_BUILD_DIR}/${MY_MODULE}/CMakeCache.txt

# THE END
ENV DEBIAN_FRONTEND teletype
