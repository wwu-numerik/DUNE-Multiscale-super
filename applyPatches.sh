#!/bin/bash

echo "Apply all patches needed for this state to compile and run (only the msfem, though!)..."
cd dune-multiscale
patch -p0 < ../patches/multiscale.patch
cd ..
echo "... done!"
exit 0

