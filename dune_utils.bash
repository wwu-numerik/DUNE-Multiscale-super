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

CMD=""
if [ $(type -t sd) ] ; then
  CMD="sd "
fi

if [ $(ionice -c 3 /bin/true) ] ; then
  CMD="${CMD} ionice -c 3  "
fi

CMD="${CMD} nice time ./dune-common/bin/dunecontrol"
