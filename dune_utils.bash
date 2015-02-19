#assumes global $OPTS
function getOptsFile( )
{
if [ x${1} = x ] ; then
  OPTS=config.opts.wwu_no_documentation
  ln -sf ${OPTS} config.opts.last
else
	OPTS=${1}
	ln -sf ${OPTS} config.opts.last
fi
}


if [ $(type -t sd) ] ; then
  CMD="sd ionice -c 3 nice time ./dune-common/bin/dunecontrol"
else
  CMD="ionice -c 3 nice time ./dune-common/bin/dunecontrol"
fi