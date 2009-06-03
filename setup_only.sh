#!/bin/bash

#!/bin/bash

if [ $# == 0 ]; then
	echo "needs module name as parameter"
else
	if [ $# == 1 ]; then
		./dune-common/bin/dunecontrol --only=$1 --opts=config.opts all
	else
		./dune-common/bin/dunecontrol --only=$1 --opts=$2 all
	fi
fi
