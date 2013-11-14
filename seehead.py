#!/usr/bin/env python
#Display information about a .fits file.
import sys

try:
	import pyfits
except:
	print 'You need PyFits to do that.'

try:
	inputFits = sys.argv[1]
except:
	print ''

#open file:
try:
	hdulist = pyfits.open(inputFits)
except:
	sys.exit('Could not read file: '+inputFits)

try:
	#requested extension:
	ext=int(sys.argv[2])
	print hdulist[ext].header
except:
	#if no requested extension, or if not present:
	hdulist.info()

hdulist.close()
