#!/usr/bin/env python
#Tristan Cantat-Gaudin, 29/07/2013.
#Converts .fits file to ascii.
#Write two columns: wavelength and flux.
#Syntax:
#	./fits2txt.py file.fits
#Output:
#	file.txt
#
#	Options:
#		-noH 	: print no header
#		-v	: verbose
#		min=, max=	: to convert only a part of the spectrum
#
#
import pyfits
import numpy as np
import sys
import array

########################### DATA ################################
try:
	fileName = sys.argv[1]
except IndexError:
	print 'You need to give a file name.'
	sys.exit()


outFile=fileName.split('.fits')[0]+'.txt'


try:
##########             Read the fits file:
	hdulist = pyfits.open(fileName)
	#hdulist.info()			#displays info about the content of the file
					#(what we use for Daospec has only ONE extension)
	#print hdulist[0].header	#to print the whole header!
	wave_base = hdulist[0].header['CRVAL1']	# Angstrom
	try:
		wave_step = hdulist[0].header['CD1_1']	# Angstrom
	except:
		wave_step = hdulist[0].header['CDELT1']	# Angstrom
	flux = hdulist[0].data
	waveobs = np.arange(wave_base, wave_base+len(flux)*wave_step, wave_step)
	if len(waveobs) == len(flux) + 1:
		waveobs = waveobs[:-1]
	hdulist.close()
except:
	print 'Problem reading the file',fileName
	sys.exit()


#read min and max wavelength from arguments:
minWav=min(waveobs)
maxWav=max(waveobs)
for opt in sys.argv:
	if 'min=' in opt:
		minWav=float(opt[4:])
	if 'max=' in opt:
		maxWav=float(opt[4:])



#and write to file:
theFile=open(outFile,"w")
if ('-noH' in sys.argv)==False:
	theFile.write('#waveobs flux\n')
for i,w in enumerate(waveobs):
	if (w>=minWav) and (w<=maxWav):
		theFile.write('%f %f %s' % (w,flux[i],'\n'))
#		theFile.write(`w`+' '+`flux[i]`+' '+`flux[i]`+' '+`flux[i]`+'\n')
theFile.close()

if '-v' in sys.argv:
	print '\twavelength step:',wave_step
	print '\tmin. wavelength:',minWav
	print '\tmax. wavelength:',maxWav
	print 'Done!'
