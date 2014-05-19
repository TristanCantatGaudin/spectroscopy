#!/usr/bin/env python
#Tristan Cantat-Gaudin, 19/05/2014.
#Degrades a spectrum from a given resolution to another.
#can read ASCII or FITS (an writes the output in the same format).
#Syntax:
#	./degrade.py file.txt 80000 47000
#Output:
#	file_deg.txt
#
# This version uses a Fortran subroutine that has to be compiled with f2py.
# If it doesn't exist, the subroutine is written out and compiled automatically 
# by the script. (f2py is distributed with NumPy)
# In some cases, "f2py -c degrade_function.f -m degradeF" won't work and
# "sudo f2py -c degrade_function.f -m degradeF" is required instead...
#
noPyFits=False
try:
	import pyfits
except:
	noPyFits=True
import numpy as np
from math import pi
from math import ceil
import sys
import os
import matplotlib.pyplot as plt

verbose=False

try:
	fileName=sys.argv[1]
except:
	print 'Syntax:'
	print './degrade.py file.txt 80000 47000'
	sys.exit()




############# READ  ###########################################
if fileName[-5:]=='.fits' or fileName[-8:-3]=='.fits':
	autoFormat='fits'
	try:	
		toto=fileName[-3:]
		if toto[0]=='[' and toto[-1]==']':
			ext=int(toto[1])
			fileName=fileName[:-3]
		else:
			ext=0
	except:
		ext=0

##########             Read the fits file:
	hdulist = pyfits.open(fileName)
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
else:
	linesToRemove=1
	#read the wavelength and flux from input txt:
	arr=[]
	try:
		inp = open(fileName,"r")
		lines = inp.readlines()[linesToRemove:]
		for line in lines:
			numbers = line.split()
			arr.append(numbers)
		waveobs = [float(el[0]) for el in arr]
		wave_step1 = (waveobs[-1]-waveobs[0])/(len(waveobs)-1)
		flux = [float(el[1]) for el in arr]
		autoFormat='ascii'
	except:
		if noPyFits:
			print 'PyFITS could not be imported. Are you trying to read a fits file?'
		sys.exit('Problem reading the spectrum \"'+fileName+'\".')














#requested change of resolution:
Ri=float(sys.argv[2])
Rf=float(sys.argv[3])

k=2.354820045

#Know how many pixels to use in the width of the gaussian (three sigmas is best):
fwi=1.*max(waveobs)/(2*Ri)
fwf=1.*max(waveobs)/(2*Rf)
sigma=((fwf**2-fwi**2)/k)**(0.5)
step = (max(waveobs)-min(waveobs))/len(waveobs)
widthGaussian=3*int(ceil( sigma/step ))
if verbose==True:
	print 'Width of Gaussian:',widthGaussian,'pixels.'


#---------------------- THIS BLOCK WILL USE A FORTRAN SUBROUTINE
dim = len(waveobs)
degradedFlux = np.zeros(dim) #set an empty array that will receive additional fluxes
totalTransferedFlux = np.zeros(dim) #for test...

try:
	from degradeF import convolve
except:
	print 'CONVOLVE subroutine not found. Compiling it on the spot.'
	fstring = 'C ----------------------------------\n      SUBROUTINE CONVOLVE(waveobs,flux,dF,wiG,Ri,Rf,N)\nC\nC    CONVOLVES THE FLUX WITH A SLIDING GAUSSIAN\nC            \nC\n      INTEGER N\n      REAL*8 Ri\n      REAL*8 Rf\n      INTEGER wiG\n      REAL*8 dF(N)\n      REAL*8 flux(N)\n      REAL*8 waveobs(N)\nC\n      REAL*8 k\n      REAL*8 pi\n      REAL*8 fwi\n      REAL*8 fwf\n      REAL*8 sigmag_c\n      REAL*8 g\n      INTEGER boundmin\n      INTEGER boundmax\n      REAL*8 ww\n      INTEGER egg\n\n      k=2.35482005\n      pi=3.1415927\n\n      DO I=1,N\n         fwi=1.*waveobs(I)/(2*Ri)\n         fwf=1.*waveobs(I)/(2*Rf)\n         sigmag_c=(fwf**2-fwi**2)/k\n         A=1./(2*pi*sigmag_c)**(0.5)\nC        boundaries for the window:\n         IF (I.GT.wiG) THEN\n            boundmin=I-wiG\n         ELSE\n            boundmin=0\n         ENDIF\n         egg=N-wiG\n         IF (I.LT.egg) THEN\n            boundmax=I+wiG\n         ELSE\n            boundmax=N-I\n         ENDIF\nC        loop inside that window:\n         DO J=boundmin,boundmax\n            ww = waveobs(J)\n            g = A*EXP( -((waveobs(I)-ww)**2)/(2*sigmag_c))\n            dF(J)=dF(J)+g*flux(I)\n         ENDDO\n      ENDDO\n      END\n'
	ffile=open('degrade_function.f','w')
	ffile.write(fstring)
	ffile.close()
	os.system('f2py -c degrade_function.f -m degradeF')
	#"f2py -c degrade_function.f -m degradeF"
	try:
		from degradeF import convolve
	except:
		print 'Problem with the compilation of the Fortran subroutine "degrade_function.f" using f2py.'
		sys.exit()

convolve( waveobs, flux, degradedFlux, widthGaussian, Ri, Rf, dim)	#called the Fotran!
#---------------------- END OF "FORTRAN BLOCK"

#normalize it:
totOrig=sum( flux )
totDeg=sum( degradedFlux )
degradedFlux=[1.*a*totOrig/totDeg for a in degradedFlux]


#cut out the edges:
waveobsDeg = waveobs[widthGaussian:-1*widthGaussian]
degradedFlux = degradedFlux[widthGaussian:-1*widthGaussian]






########### WRITE THE RESULT TO A FILE #############
if autoFormat=='fits':
	outputSpectrum=fileName.replace('.fits','_deg.fits')
	#create the fits:
	flux = np.array(degradedFlux,dtype='float32') #important for Daospec to have float32!!!
	os.system('rm -f '+outputSpectrum)
	pyfits.writeto(outputSpectrum,flux)

	header = pyfits.getheader(outputSpectrum)
	header.update('CRVAL1', min(waveobsDeg), "wavelength zeropoint")
	header.update('CD1_1', wave_step, "wavelength step")
	header.update('CDELT1', wave_step, "wavelength step")
	header.update('CRPIX1', 1.0, "Pixel zeropoint")
	header.update('NAXIS', 1, "Number of axes")
	header.update('NAXIS1', len(flux), "Axis length")

	os.system('rm -f '+outputSpectrum)
	pyfits.writeto(outputSpectrum,flux,header)

elif autoFormat=='ascii':
	typ=fileName.split('.')[-1]
	outFile=fileName.replace('.'+typ,'_deg.'+typ)
	theFile=open(outFile,"w")
	theFile.write('#waveobs flux\n')
	for i,w in enumerate(waveobsDeg):
		theFile.write(`w`+' '+`degradedFlux[i]`+'\n')
	theFile.close()

if verbose==True:
	print 'Done.'
