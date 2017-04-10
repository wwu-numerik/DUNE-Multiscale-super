#/bin/bash!
echo "Usage: <total num of p> , <p fast>, <p med.>, <p slow> <gird sz> <corr len> <sigma> <levels> [out file]"
mpirun -np ${1} ./dune-random-diff/build-release/src/dune-random-diff ${2} ${3} ${4} ${5} 1 ${6} ${7} 5e-3 3 ${8} ${9}
