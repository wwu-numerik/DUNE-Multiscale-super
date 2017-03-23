#assumes global $OPTS
#prefix for the root of the project
#By default base path is assumed to be current working directory
basePath=$(dirname -- $(readlink -fn -- "$0"))

#echo some infromation
echo -e "Base path for the dune modules:$basePath\n";

echo -e  "Checking directory structure...\n"
if [ ! -d "$basePath/$varDune" ]; then
  echo -e "Could not find dune-multiscale-super\n"
fi

path="$basePath/$varDune"


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

if  [ $(ionice -c 3 /bin/true) ] ; then
  IONICE="ionice -c 3"
else
  IONICE=""
fi

if [ $(type -t sd) ] ; then
  CMD="sd ${IONICE} nice time $path/dune-common/bin/dunecontrol"
else
  CMD="${IONICE} nice $path/dune-common/bin/dunecontrol"
fi
