#!/usr/bin/env python
#Tristan Cantat-Gaudin, 14/03/2013.
#Degrades a spectrum from a given resolution to another.
#can read ASCII or FITS (an writes the output in the same format).
#Syntax:
#	./degrade.py file.txt 80000 47000
#Output:
#	file_deg.txt
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
if fileName[-5:]=='.fits':
	autoFormat='fits'
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
		if noPyfits:
			print 'PyFITS could not be imported. Are you trying to read a fits file?'
		sys.exit('Problem reading the spectrum \"'+fileName+'\".')














#requested change of resolution:
Ri=float(sys.argv[2])
Rf=float(sys.argv[3])

k=2.354820045

#Know how many pixels to use in the width of the gaussian (two sigmas):
fwi=1.*max(waveobs)/(2*Ri)
fwf=1.*max(waveobs)/(2*Rf)
sigma=((fwf**2-fwi**2)/k)**(0.5)
step = (max(waveobs)-min(waveobs))/len(waveobs)
widthGaussian=2*int(ceil( sigma/step ))
if True:
	print 'Width of Gaussian:',widthGaussian,'pixels.'



degradedFlux=[0 for el in waveobs] #set an empty array that will receive additional fluxes
totalTransferedFlux=[0 for el in waveobs] #for test...
for i,w in enumerate(waveobs):
	#make a gaussian of requested width:
	fwi=1.*w/(2*Ri)
	fwf=1.*w/(2*Rf)
	sigmag_c=(fwf**2-fwi**2)/k
	A=1./(2*pi*sigmag_c)**(0.5)
	#print 'lambda:',w,'sigma:',sigmag_c
	print round(100.*i/len(waveobs),2),'% done.'
	def gauss_0(x):
		global A,w,sigmag_c
		g = A*np.exp( -((x-w)**2)/(2*sigmag_c))
		return g


	#set the boundaries of the range in which to pick flux:
	if i>widthGaussian:
		boundmin=i-widthGaussian
	else:
		boundmin=0
	if i<len(waveobs)-widthGaussian:
		boundmax=i+widthGaussian
	else:
		boundmax=len(waveobs)-1

	if verbose:
		print 'The original resolution is',Ri
		print 'At this resolution, for a wavelength of',w
		print 'the fw to distinguish between two lines is:',fwi
		print 'The final resolution is',Rf
		print 'so the final resolved lines will have a fw of',fwf
		print 'and the gaussian to convolve with has sigma:',(sigmag_c)**(0.5)
		print 'We are at index',i,', picking flux in index range:',boundmin,'to',boundmax,'.'

	#and pick flux in that range:
	plt.figure(0)
	for j in np.arange(boundmin,boundmax+1):
		w2=waveobs[j]
		#totalTransferedFlux[i]=totalTransferedFlux[i]+gauss_0(w2)*flux[i]
		#print '   flux is:',flux[i],'at w2=',w2,'(i=',i,'gaussian gives:',gauss_0(w2)
		#print '       so bin #',j,'receives +',gauss_0(w2)*flux[i]
		degradedFlux[j]=degradedFlux[j]+gauss_0(w2)*flux[i]
	#print 'total transferred:',totalTransferedFlux[i]


#normalize it:
totOrig=sum( flux )
totDeg=sum( degradedFlux )
degradedFlux=[1.*a*totOrig/totDeg for a in degradedFlux]

#for tests ###############
#diff=[a-b for a,b in zip(degradedFlux,flux)]
#totalTransferedFlux=[1.*a/b for a,b in zip(totalTransferedFlux,flux)]
#degCorr=[a/b for a,b in zip(degradedFlux,totalTransferedFlux)]
#test=[a/b for a,b in zip(flux,degCorr)]
##########################

#cut out the edges:
waveobsDeg = waveobs[widthGaussian:-1*widthGaussian]
degradedFlux = degradedFlux[widthGaussian:-1*widthGaussian]

if False:
	plt.figure(0)
	plt.plot(waveobsDeg,degradedFlux,'g-')
	plt.plot(waveobs,flux,'b-')
	#plt.plot(waveobs,diff,'k-')
	#plt.plot(waveobs,totalTransferedFlux,'r-')
	#plt.plot(waveobs,degCorr,'y-')
	#plt.ylim(-10,10000)
	plt.show()







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
	ext=fileName.split('.')[-1]
	outFile=fileName.replace('.'+ext,'_deg.'+ext)
	theFile=open(outFile,"w")
	theFile.write('#waveobs flux\n')
	for i,w in enumerate(waveobsDeg):
		theFile.write(`w`+' '+`degradedFlux[i]`+'\n')
	theFile.close()
