#/bin/bash!

####################################################################################################
# Bash script for helping mentain and build dune-multiscalesuper set of libs on ITWM beehivecluster
# Exports paths for dune-multiscalesuper block
# verssion: 2.0
####################################################################################################

#prefix for the root of the project
#By default base path is assumed to be current working directory
basePath=$(dirname -- $(readlink -fn -- "$0"))
varOpt="./opt"
varDune="."

if [ "$#" -eq 1 ]; then
  echo "Setting dir name to user defined: $1"
  basePath=$1
fi
#echo some infromation
echo -e "Base path for the dune modules:$basePath\n";

echo -e  "Checking directory structure...\n"
if [ ! -d "$basePath/$varOpt" ] || [ ! -d "$basePath/$varDune" ]; then
  echo -e "Could not find dune or/and opt\n"
  exit 1
fi

echo -e "Every external dep. is to be put in a separte directory under ./opt. If new dep. are needed it should be putted in ./opt in separate directory and noted in the script\n"

LIBS=""
INCLUDES=""


#OPEN BLAS
####################################################################################################
echo -e "Exporting openblas library"
varBLAS="$basePath/$varOpt/blas"
if [ ! -d $varBLAS ]; then
  echo "Could not find OpenBlas dir"
fi
PATH=$varBLAS/bin:${PATH}
LIBS=$varBLAS/lib:${LIBS}
INCLUDES=$varBLAS/include:${INCLUDES}
####################################################################################################

#SuiteSparce
####################################################################################################
echo -e "Exporting SuiteSparce library"
varSuiteSparce="$basePath/$varOpt/SS4"
if [ ! -d $varSuiteSparce ]; then
  echo "Could not find SuiteSparce dir"
fi

PATH=$varSuiteSparce/bin:${PATH}
LIBS=$varSuiteSparce/lib:${LIBS}
INCLUDES=$varSuiteSparce/include:${INCLUDES}
####################################################################################################


#SuperLu_ROOT
####################################################################################################
echo -e "Exporting SuperLU library root"
varSuperLu="$basePath/$varOpt/SuperLU_4.3"
if [ ! -d $varSuperLu ]; then
  echo "Could not find SuperLU dir"
fi
LIBS=$varSuperLu/lib:${LIBS}
INCLUDES=$varSuperLu/SRC:${INCLUDES}

####################################################################################################


#FFTW library
####################################################################################################
echo -e "Exporting libfftw3"
varFFT="$basePath/$varOpt/libfftw3"
if [ ! -d $varFFT ]; then
  echo "Could not find libfftw3 library dir"
fi

PATH=$varFFT/bin:${PATH}
LIBS=$varFFT/lib:${LIBS}
INCLUDES=$varFFT/include:${INCLUDES}
####################################################################################################

#ParMetis
####################################################################################################
echo -e "Exporting ParMetis Lib library"
varParMetis="$basePath/$varOpt/ParMetis4"
if [ ! -d $varTbb ]; then
  echo "Could not find ParMetisLib library dir"
fi

PATH=$varParMetis/bin:${PATH}
LIBS=$varParMetis/lib:${LIBS}
INCLUDES=$varParMetis/include:${INCLUDES}
####################################################################################################

#TBB library from intel
####################################################################################################
echo -e "Exporting intel tbb library"
varTbb="$basePath/$varOpt/tbb"
if [ ! -d $varTbb ]; then
  echo "Could not find TBB library dir"
fi

PATH=$varTbb/bin:${PATH}
LIBS=$varTbb/lib:${LIBS}
INCLUDES=$varTbb/include:${INCLUDES}
####################################################################################################

#export PATH
#export LD_LIBRARY_PATH=${LIBS}:${LD_LIBRARY_PATH}
#export C_INCLUDE_PATH=${INCLUDES}:${C_INCLUDE_PATH}
#export CPLUS_INCLUDE_PATH=${INCLUDES}:${CPLUS_INCLUDE_PATH}

#Python 
####################################################################################################
#varPy="$basePath/$varOpt/py"
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
module add lang/python/3.5.0

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

echo -e "\nExporting Dune Contorol path envirment variable\n"
export DUNE_CONTROL_PATH=${basePath}/${varDune}

echo -e "\nExporting Dune dep paths for cmake...\n"

# Export real paths
##################################################################################
export DUNE_BLAS=$(cd ${varBLAS}; pwd)
export DUNE_FFTW_ROOT=$(cd ${varFFT}; pwd)
export DUNE_SUITESPARSE_LIBRARY_ROOT=$(cd ${varSuiteSparce}; pwd) 

export DUNE_SUPERLU_ROOT=$(cd ${varSuperLu}; pwd) 
export DUNE_PARMETIS_ROOT=$(cd ${varParMetis}; pwd) 

export DUNE_TBB_ROOT=$(cd ${varTbb}; pwd) 


# Export real paths to Boost
#################################
#export ITWM_BOOST_LINKER_PATH="/m/soft/boost/boost_1_61_0/lib"
export DUNE_BOOST_ROOT="/m/soft/boost/boost_1_61_0"

##################################################################################

# Should not be needed (cmake gives r-path)
##################################################################################
#export PATH
export LD_LIBRARY_PATH=${LIBS}:${LD_LIBRARY_PATH}
#export C_INCLUDE_PATH=${INCLUDES}:${C_INCLUDE_PATH}
#export CPLUS_INCLUDE_PATH=${INCLUDES}:${CPLUS_INCLUDE_PATH}

echo -e "Done..."



