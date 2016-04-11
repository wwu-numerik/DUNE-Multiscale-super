#!/bin/bash

set -e

THISDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
BASEDIR=${TRAVIS_BUILD_DIR:-${THISDIR}/..}
ID=${TRAVIS_TAG:-${TRAVIS_BRANCH}}

pushd ${BASEDIR}

git clone git@github.com:wwu-numerik/wwu-numerik.github.io.git site

BUILDDIR=${PWD}/build
${THISDIR}/build_docs.sh ${BUILDDIR}

cd site
git config user.name "DUNE Community Bot"
git config user.email "dune-community.bot@wwu.de"

for i in multiscale mlmc; do
  TARGET=docs/dune-${i}/${ID}/
  mkdir -p ${TARGET}
  rsync -a --delete  ${BUILDDIR}/dune-${i}/doc/doxygen/html/ ${TARGET}/
  git add ${TARGET}
done

git commit -m "Updated documentation for dune-multiscale/mlmc ${ID}"
git push

popd