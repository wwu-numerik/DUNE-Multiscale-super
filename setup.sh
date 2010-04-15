#!/bin/bash

if [ x$1 = x ] ; then
	OPTS=config.opts.wwu_no_documentation
else
	OPTS=${1}
fi
ln -sf ${OPTS} config.opts.last
./dune-common/bin/dunecontrol --opts=${1} all