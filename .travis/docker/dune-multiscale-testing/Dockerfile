FROM dunecommunity/testing-base:multiscale_pre_2.4

# COMPILER Selection (ARGS cannot be uppercase)
ARG cc
ARG cxx
ENV CC $cc
ENV CXX $cxx
ENV MY_MODULE dune-multiscale
ENV OPTS config.opts/travis.ninja
ENV PATH /usr/lib/ccache:${PATH}
ENV SRC_DCTRL ${SUPERDIR}/dune-common/bin/dunecontrol

WORKDIR ${SUPERDIR}
ADD . ${SUPERDIR}

RUN cd ${SUPERDIR} && \
    git submodule foreach git reset --hard && \
    git submodule update --init --recursive && \
    ./local/bin/download_external_libraries.py && \
    ./local/bin/build_external_libraries.py && \
    ${SRC_DCTRL} ${BLD} all && \
    rm -rf ${MY_MODULE} && mkdir ${MY_MODULE}
