#!/bin/bash

echo "Apply all patches needed for this state to compile and run (only the msfem, though!)..."
cd dune-spgrid
patch -p0 < ../patches/spgrid.patch
cd ../dune-stuff
patch -p0 < ../patches/stuff.patch
cd ../dune-subgrid
patch -p0 < ../patches/subgrid.patch
cd ../dune-multiscale
patch -p0 < ../patches/multiscale.patch
cd ..
echo "... done!"
exit 0

