#!/bin/bash
set -e
. dune_utils.bash
getOptsFile $1
${CMD} --opts=${1} all
