#!/usr/bin/env bash
set -e
cd $(dirname $0)


BN=$(basename ${PWD})
REV=$(git describe --always --dirty)
TARFILE=${BN}-${REV}.tar
rm -f ${TARFILE}

./git-archive-all --prefix ${BN}-${REV}/  ${TARFILE}
cat ${TARFILE} | gzip > ${TARFILE}.gz
cat ${TARFILE} | bzip2 > ${TARFILE}.bz2
rm ${TARFILE}
