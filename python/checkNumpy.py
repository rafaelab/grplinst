# Retrieve NumPy include path.
# If numpy can't be imported, return -1 and exit without throwing an exception.

import sys

try:
	import numpy
	from pkg_resources import parse_version
	if parse_version(numpy.__version__) < parse_version('1.6.0'):
		sys.exit(-1)
	sys.stdout.write(numpy.get_include())
except ImportError:
	sys.exit(-1)