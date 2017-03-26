#/bin/bash!

####################################################################################################
# Bash script for helping mentain and build dune-multiscale-super set of libs on ITWM beehivecluster
# Exports paths for dune-multiscale-super block
# verssion: 2.1
# TO DO: ADD LIBS AS A FUNCTION
####################################################################################################

#prefix for the root of the project
#By default base path is assumed to be current working directory
basePath=$(dirname -- $(readlink -fn -- "$0"))

if [ "$#" -eq 1 ]; then
  echo "Setting dir name to user defined: $1"
  basePath=$(cd ${1}; pwd)
fi
#echo some infromation
echo -e "Base path for optional libs:${baseOpt}\n"
echo -e "Checking directory structure...\n"
if [ ! -d ${basePath} ]; then
  echo -e "Could not find specified dir"
  exit 1
fi

echo -e "Every external dep. is to be put in a separte directory under ${baseOpt} If new dep. are needed it should be putted in ${baseOpt} in separate directory and noted in this script\n"


#Open Blas
####################################################################################################
echo -e "Exporting openblas library"
varBLAS="$basePath/blas"
if [ ! -d $varBLAS ]; then
  echo "Could not find OpenBlas dir"
fi
####################################################################################################


#SuiteSparce
####################################################################################################
echo -e "Exporting SuiteSparce library"
varSuiteSparce="$basePath/SS4"
if [ ! -d $varSuiteSparce ]; then
  echo "Could not find SuiteSparce dir"
fi
####################################################################################################


#SuperLU
####################################################################################################
echo -e "Exporting SuperLU library root"
varSuperLu="$basePath/SuperLU_4.3"
if [ ! -d $varSuperLu ]; then
  echo "Could not find SuperLU dir"
fi

####################################################################################################


#FFTW library
####################################################################################################
echo -e "Exporting libfftw3"
varFFT="$basePath/libfftw3"
if [ ! -d $varFFT ]; then
  echo "Could not find libfftw3 dir"
fi
####################################################################################################


#ParMetis
####################################################################################################
echo -e "Exporting ParMetis Lib library"
varParMetis="$basePath/ParMetis4"
if [ ! -d $varTbb ]; then
  echo "Could not find ParMetis dir"
fi
####################################################################################################


#TBB
####################################################################################################
echo -e "Exporting intel tbb library"
varTbb="$basePath/tbb"
if [ ! -d $varTbb ]; then
  echo "Could not find TBB dir"
fi
####################################################################################################

#Python 
####################################################################################################
#varPy="$basePath/py"
#if [ ! -d $varPy ]; then
#  echo "Could not find python dir"
#fi
#
#echo -e "Searching for: lib/python3.5/site-packages"
#if [ ! -d $varPy/lib/python3.5/site-packages ]; then
#  echo "Could not find python site-packages"
#fi
#
#echo -e "Exporting pyhton libs"
#
#export PYTHONPATH=$varPy/lib/python3.5/site-packages:${PYTHONPATH}
#export PATH=$varPy/bin:${PATH}
####################################################################################################

echo -e "Setting up the modules...\n"
#Additional libs
module add mpi/openmpi-x86_64
module add lib/boost/1_61_0
module add lang/python/3.5.2

#IDEs and compiler tools
module add develop/qtcreator
module add develop/qt
module add devel/cmake/3.6.1

#Addtitional software
module add soft/paraview
module add soft/texlive

#Compilers
module add compiler/gcc/5.4.0

#Print loaded modules on the system
module list

echo -e "\nExporting Dune dep paths for cmake...\n"

function exportDune
{
  if [ -d ${1} ]; then
    export ${2}=$(cd ${1}; pwd)
  fi
}

# Export real paths
##################################################################################

exportDune ${varBLAS} DUNE_BLAS
exportDune ${varSuiteSparce} DUNE_SUITESPARSE_ROOT
exportDune ${varSuperLu} DUNE_SUPERLU_ROOT
exportDune ${varParMetis} DUNE_PARMETIS_ROOT
exportDune ${varFFT} DUNE_FFTW_ROOT
exportDune ${varTbb} DUNE_TBB_ROOT

# Export real paths to Boost
#################################
#export ITWM_BOOST_LINKER_PATH="/m/soft/boost/boost_1_61_0/lib"
export DUNE_BOOST_ROOT="/m/soft/boost/boost_1_61_0"
##################################################################################

echo -e "Done..."
