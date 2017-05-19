DIR="$(cd "$(dirname ${BASH_SOURCE[0]})" ; cd ../../ ; pwd -P )"

docker run -it -v ${DIR}/dune-multiscale:/root/src/dune-multiscale \
    -v ${DIR}/dune-mlmc:/root/src/dune-mlmc \
    dunecommunity/dune-multiscale-testing_gcc:master
