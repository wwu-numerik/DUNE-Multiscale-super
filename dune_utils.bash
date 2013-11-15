#assumes global $OPTS
function getOptsFile( )
{
if [ x${1} = x ] ; then
	if [ -e config.opts.last ] ; then
		OPTS=$(readlink config.opts.last)
	else
		OPTS=config.opts.wwu_no_documentation
		ln -sf ${OPTS} config.opts.last
	fi
else
	OPTS=${1}
	ln -sf ${OPTS} config.opts.last
fi
}
#export LD_LIBRARY_PATH=/Applic.HPC/intel/impi/4.1.0.024/lib64:/Applic.HPC/boost/boost_1_53_0/stage/lib:/Applic.HPC/gcc/4.7.2/lib64:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/Applic.HPC/openmpi/gcc_mellanox/1.6.5/lib/:/Applic.HPC/boost/boost_1_53_0/stage/lib:/Applic.HPC/gcc/4.7.2/lib64:$LD_LIBRARY_PATH
#export BOOST_LIB_DIR=/Applic.HPC/boost/boost_1_53_0/stage/lib
#export BOOST_ROOT=/Applic.HPC/boost/boost_1_53_0
#module load compiler/gcc/4.7.2 tools/mic lib/boost/1.53 compiler/intel/14 mpi/intel/4.1.0.024
module load compiler/gcc/4.7.2 tools/mic lib/boost/1.53 compiler/intel/14 mpi/openmpi/gcc/1.6.5


