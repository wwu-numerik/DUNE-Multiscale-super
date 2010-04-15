#!/bin/bash

branch=$(git branch --no-color | cut -f 2 -d " ")

for postfix in grid grid-howto common istl; do
	i=dune-$postfix
	cd ${i}
	git fetch origin
	new=$(git rev-parse origin/$branch)
	git checkout $new
	cd ..
	git add $i
done
