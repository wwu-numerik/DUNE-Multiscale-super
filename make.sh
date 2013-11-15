#!/bin/bash
set -e
. dune_utils.bash
getOptsFile $2
./dune-common/bin/dunecontrol --module=$1 --opts=${OPTS} --use-cmake make -k
