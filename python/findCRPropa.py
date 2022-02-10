import sys

try:
	import crpropa
except ImportError:
	sys.exit(-1)

if sys.argv[1] == 'swig_interface':
	sys.stdout.write(crpropa.getDataPath('swig_interface'))
elif sys.argv[1] == 'install_prefix':
	sys.stdout.write(crpropa.getInstallPrefix())
