#!/bin/bash

if [ x$2 = x ] ; then
	OPTS=config.opts.wwu_no_documentation
else
	OPTS=${2}
fi
ln -sf ${OPTS} config.opts.last
./dune-common/bin/dunecontrol --only=$1 --opts=${OPTS} all
