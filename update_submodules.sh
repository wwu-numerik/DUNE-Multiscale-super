#!/bin/bash

branch=$(git branch --no-color | cut -f 2 -d " ")

for i in $(git submodule status | cut -f 2 -d " " ) ; do
	cd ${i}
	git fetch origin
	new=$(git rev-parse origin/$branch)
	git checkout $new
	cd ..
	git add $i
done

git commit -a -m "submodules updated to latest"
git push origin $branch

