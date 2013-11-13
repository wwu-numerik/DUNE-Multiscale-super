#PBS -o dunecontroloutput.dat
#PBS -l walltime=00:15:00,nodes=1:ppn=8
#PBS -A o0num
# #PBS -M s_kaul01@uni-muenster.de
#PBS -m ae
#PBS -q mescashort
#PBS -N dunecontrol
#PBS -j oe
cd /scratch/tmp/s_kaul01/dune-ms-super
export LD_LIBRARY_PATH=/Applic.PALMA/compiler/gcc/4.7.2/lib64:$LD_LIBRARY_PATH
# export BOOST_ROOT=/Applic.PALMA/boost/boost_1_48_0
# source BOOST_PATH.sh

# this should build dune-common dune-geometry dune-grid dune-istl dune-localfunctions dune-fem dune-subgrid dune-stuff dune-spgrid dune-multiscale
./dune-common/bin/dunecontrol --use-cmake --opts=config.opts.palma --only=dune-fem all
