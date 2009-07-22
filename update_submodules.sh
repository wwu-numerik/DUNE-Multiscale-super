#!/bin/bash

for i in $(find . -mindepth 1  -maxdepth 1 -type d ) ; do 
	cd ${i}
	git fetch origin
	git merge origin/master
	cd ..
done

git commit -a -m"submodules updated to latest"
git push origin master

