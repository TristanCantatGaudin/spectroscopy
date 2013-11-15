#!/usr/bin/env python
#cut a slice of a spectrum:
#
#    ./cutSpec.py spec.fits 4800 5800 specL.fits
#    ./cutSpec.py spec.fits[1] 4800 5800 specL.fits
#
#NB: the output will be in 32 bit, whatever the input!
#    the code can easily be modified to output 64bits.


try:
	import pyfits
	noPyfits=False
except:
	noPyfits=True
	sys.exit('You need PyFits to do that.')
import numpy as np
import sys

overwrite=True	#this sets clobber=True when pyfits.writeto writes the final file





# Read arguments:
try:
	infile = sys.argv[1]
	minWAV = float(sys.argv[2])
	maxWAV = float(sys.argv[3])
	outfile = sys.argv[4]
except IndexError:
	print 'The syntax is:'
	print sys.argv[0], "file.txt -options"
	sys.exit()

#find in requested a specific extension, ext=0 if not:
try:	
	ext=infile[-3:]
	if ext[0]=='[' and ext[-1]==']':
		ext=int(ext[1])
	else:
		ext=0
except:
	ext=0





#Read file:
try:
	hdulist = pyfits.open(infile)
except:
	sys.exit('Could not read file: '+infile)
#hdulist.info()			#displays info about the content of the file
				#(what we use for Daospec has only ONE extension)
#print hdulist[0].header	#to print the whole header!
wave_base = hdulist[ext].header['CRVAL1']	# Angstrom
try:
	wave_step = hdulist[ext].header['CD1_1']	# Angstrom
except:
	wave_step = hdulist[ext].header['CDELT1']	# Angstrom
flux = hdulist[ext].data
waveobs = np.arange(wave_base, wave_base+len(flux)*wave_step, wave_step)
if len(waveobs) == len(flux) + 1:
	waveobs = waveobs[:-1]
hdulist.close()


#Now keep only the requested wavelength range:
cutFlux = [ f for f,w in zip(flux,waveobs) if w>minWAV and w<maxWAV ]
cutFlux = np.array(cutFlux,dtype='float32')	# <---- sets the data to 32bits format

cutWaveobs = [ w for w in waveobs if w>minWAV and w<maxWAV ]


#update the header and write the fits file:
theHeader = pyfits.getheader(infile)
theHeader.update('CRVAL1', min(cutWaveobs), "wavelength zeropoint")
theHeader.update('NAXIS1', len(cutFlux), "Axis length")
theHeader.add_history('')
theHeader.add_history(' Reduced the spectral range:')
theHeader.add_history(' from ['+str(min(waveobs))+'-'+str(max(waveobs))+']')
theHeader.add_history(' to ['+str(minWAV)+'-'+str(maxWAV)+'].')
pyfits.writeto(outfile,cutFlux,theHeader,clobber=overwrite)
