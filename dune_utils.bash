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
