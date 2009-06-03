#!/bin/bash

if [ $# == 0 ]; then
	./dune-common/bin/dunecontrol --opts=config.opts all
else
	./dune-common/bin/dunecontrol --opts=$1 all	
fi
