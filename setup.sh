#!/bin/bash
set -e
. dune_utils.bash
getOptsFile $1
./dune-common/bin/dunecontrol --opts=${1} --use-cmake all
