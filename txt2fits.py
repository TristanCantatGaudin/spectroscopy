#!/usr/bin/env python
#Tristan Cantat-Gaudin, 10/06/2013
#Conversion of an ASCII table
#into a FITS spectrum.
#syntax: ./txt2fits file.txt -options
#
#   options: 	-h2	to ignore 2 lines of header
#		-m10	to multiply the wavelength by 10 (if you had nm and want Angstroms in the FITS)
#		-gaps	to fill possible gaps with random noise
#		-spikes	to smooth out spikes
#		-n	to avoid negative fluxes
#		maxflux=1.2	to cut flux above 1.2
#		-v	to make the script talkative
#		ext=_deg.fits	to rename Foo.txt in Foo_deg.fits
#		cdelt=0.0147	to resample with a step of 0.0147*
#		rv=150	to add a 150km/h radial velocity (this increases the wavelength)
#		-col3	to read flux from column 3 instead of 1 (nb: Python indices start at 0!)
#
#		*fits spectra assume a constant wavelength step, so if your original ascii is not on a
#		 constant step you NEED to resample.
#
#NB: this assumes that the first column of the ASCII file is the wavelength.
#
import sys
import os
import pyfits
import numpy
from scipy.interpolate import interp1d

try:
	infile = sys.argv[1]
	
except IndexError:
	print 'The syntax is:'
	print sys.argv[0], "file.txt -options"
	sys.exit()

#read options:
linesToRemove=0
multiFactor=1.0
fillGaps=False
fillNegatives=False
maxFlux=False
cutSpikes=False
verbose=False
options = sys.argv[2:]
extension='.fits'
resample=False
radvel=0
colFlux=1
for el in options:
	if '-h' in el:
		linesToRemove = int(el.split('-h')[1])
	if '-m' in el:
		multiFactor=float(el.split('-m')[1])
	if '-gaps' in el:
		fillGaps=True
	if '-spikes' in el:
		cutSpikes=True
		try:
			threshold=float(el.split('-spikes')[1])
		except:
			threshold=3
	if '-v' in el:
		verbose=True
	if '-n' in el:
		fillNegatives=True
	if 'maxFlux=' in el:
		maxFlux=float(el[8:])
	if 'ext=' in el:
		extension=el[4:]
	if 'cdelt=' in el:
		resample=True
		wave_step=float(el[6:])
	if 'rv=' in el:
		radvel=float(el[3:])
	if '-col' in el:
		colFlux=int(el[4:])
		



#read the wavelength and flux from input txt:
arr=[]
try:
	inp = open(infile,"r")
	lines = inp.readlines()[linesToRemove:]
	if (linesToRemove==1) and (verbose==True):
		print 'Ignoring the first line.'
	if (linesToRemove>1) and (verbose==True):
		print 'Ignoring the first',linesToRemove,'lines.'
	for line in lines:
		numbers = line.split()
		arr.append(numbers)
	waveobs = [float(el[0]) for el in arr]
	flux = [float(el[colFlux]) for el in arr]
except:
	sys.exit('Problem reading the ascii spectrum \"'+infile+'\".')

#################################################################
############# RESAMPLE BY INTERPOLATION IF ASKED !!! ############
if resample==True:
	if verbose==True:
		print 'Interpolating...'
	fluxNew = interp1d(waveobs,flux,kind='linear')
	flux_new=[]
	wave_new = [waveobs[0]+k*wave_step for k in range(0,int((waveobs[-1]-waveobs[0])/wave_step))]
	for w_new in wave_new:
		#print 'min',waveobs[0]
		#print 'max',waveobs[-1]
		flux_new.append( fluxNew(w_new) )
		#print w_new, fluxNew(w_new)
	waveobs=wave_new
	flux=flux_new
#################################################################
#################################################################

#multiply (in case requested)
waveobs = [el*multiFactor for el in waveobs]
if (multiFactor!=1) and (verbose==True):
	print 'Multiplying the wavelength by a factor:',multiFactor

#add some radial velocity:
if radvel!=0:
	factor_temp=1+(radvel/299792.458)
	waveobs = [el*factor_temp for el in waveobs]
	if verbose==True:
		print 'Added a radial velocity of',radvel,'km/h.'

#output name:
outputSpectrum=infile.split('.txt')[0]+extension


#calculate CDELT1 and CRVAL1:
wave_base = waveobs[0]
if resample==False:
	wave_step = (waveobs[-1]-waveobs[0])/(len(waveobs)-1)


#fill the gaps!
if fillGaps==True:
	countGaps=0
	import random
	for i,y in enumerate(waveobs):
		if i>0:
			if waveobs[i]-waveobs[i-1] > 2*wave_step:
				countGaps=countGaps+1
				if verbose==True:
					print 'Filling a gap between '+str(waveobs[i-1])+' and '+str(waveobs[i])+'.'
				#keep the list up to flux[i-1] (included)
				list1=flux[:i]
				#make a list starting from flux[i] included
				list2=flux[i:]
				#append as many elements as needed to the first list
				nbToAdd = int((waveobs[i]-waveobs[i-1])/wave_step - 1)
				listGap=[flux[i-1]+(flux[i]-flux[i-1])*random.random() for el in range(nbToAdd)]
				#append list2 to list1
				flux=list1+listGap+list2
	wave_step = (waveobs[-1]-waveobs[0])/(len(flux)-1)
	if verbose==True:	
		print countGaps,'gap(s) filled.'


#limit max flux!
print 'Maxflux:',maxFlux
if maxFlux!=False:
	maxFluxCounter=0
	for i,el in enumerate(flux):
		if (el>maxFlux):
			flux[i]=maxFlux
			maxFluxCounter=maxFluxCounter+1
	if verbose==True:
		print maxFluxCounter,'pixels with flux higher than',maxFlux,'were found and replaced.'


#remove negative fluxes!
if fillNegatives==True:
	import random
	negativeCounter=0
	median=numpy.median(flux)
	sigma=numpy.std(flux)
	for i,el in enumerate(flux):
		if (el<0) or (el==0):
			essai=median+0.3*sigma*(random.random()-0.5)
			if essai>0:
				flux[i]=essai
			else:
				flux[i]=0.1
			negativeCounter=negativeCounter+1
	if verbose==True:
		print negativeCounter,'negative pixels found and replaced.'



#cut the spikes!
if cutSpikes==True:
	widthOfSlice=300
	spikeCounter=0
	import random
	nbSlices=int(len(flux)/widthOfSlice)
	for n in range(nbSlices):
		thisSlice=flux[widthOfSlice*n:widthOfSlice*(n+1)]
		median=numpy.median(thisSlice)
		sigma=numpy.std(thisSlice)
		if max(thisSlice) > median+threshold*sigma:
			spikeCounter=spikeCounter+1
		for i,el in enumerate(thisSlice):
			if el > median+threshold*sigma:
				flux[widthOfSlice*n+i]=median+0.3*sigma*(random.random()-0.5)
	thisSlice=flux[widthOfSlice*nbSlices:]
	median=numpy.median(thisSlice)
	sigma=numpy.std(thisSlice)
	if max(thisSlice) > median+threshold*sigma:
		spikeCounter=spikeCounter+1
	for i,el in enumerate(thisSlice):
		if el > median+threshold*sigma:
			flux[widthOfSlice*nbSlices+i]=median+0.3*sigma*(random.random()-0.5)
	if verbose==True:
		print spikeCounter,'spike(s) removed, using threshold of',threshold,'sigmas.'


#create the fits:
flux = numpy.array(flux,dtype='float32') #important for Daospec to have float32!!!
os.system('rm -f '+outputSpectrum)
pyfits.writeto(outputSpectrum,flux)

header = pyfits.getheader(outputSpectrum)
header.update('CRVAL1', wave_base, "wavelength zeropoint")
header.update('CD1_1', wave_step, "wavelength step")
header.update('CDELT1', wave_step, "wavelength step")
header.update('CRPIX1', 1.0, "Pixel zeropoint")
header.update('NAXIS', 1, "Number of axes")
header.update('NAXIS1', len(flux), "Axis length")

os.system('rm -f '+outputSpectrum)
pyfits.writeto(outputSpectrum,flux,header)
