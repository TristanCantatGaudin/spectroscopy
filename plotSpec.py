#!/usr/bin/env python
#Tristan Cantat-Gaudin, 25/07/2013.
#Quickly plot a .fits or ascii spectrum. You can also batch-plot (see below).
#
#	The available options are:
#
#		-color0=black		specify a color for spectrum 0
#					(any Matplotlib color name or hexadecimal value)
#					(default cycle: blue-cyan-grey-black-blue etc)
#
#		-rv0=-56		in km/s. Positive values will redshift the spectrum.
#					(use negative values to CORRECT for a redshift)
#
#		-label0=Sun		the label will be added in a caption.
#		-label0="Sun 1"		(use quotes if you need spaces)
#
#		-ext0=1			for fits spectra only, tells which extension to read
#					(default is 0)
#
#		-column0=1		for ascii spectra, which column contains the flux
#					(default is 1, and the wavelength has to be column 0)
#
#	/!\ REMINDER: Python indices start at 0, so the first spectrum, the first column etc
#		      are always number 0!
#
#	BATCH-PLOTTING: If you create a file "foo.bar" of instructions like:
#			star1.fits model1.txt -color0=red -color1=blue
#			star2.fits model2.txt -color0=red -color1=blue
#			star3.fits model3.txt -color0=red -color1=blue
#			You can batch-plot by calling:
#				./seefits @foo.bar
#			Type ENTER in the prompt to clear the window and execute
#			the next line of instructions.
#			Type b to go back to the previous line.
#			Type q to exit.
#
#	EXAMPLES:
#		./seefits.py file.fits
#	
#		./seefits.py file.txt
#	
#		./seefits.py file0.fits file1.txt -color0=red -color1=blue -rv1=23.7 -label1=shifted
#
try:
	import pyfits
	noPyfits=False
except:
	noPyfits=True
import matplotlib.pyplot as plt
import numpy as np
import sys
c=299792.458

#function that takes the arguments, reads the files and plots:
def actualseefits(argumentsList):
	#read the arguments:
	options=[]
	files=[]
	for arg in argumentsList[1:]:
		if arg[0]=='-':
			options.append(arg)
		else:
			files.append(arg)
	nbspec=len(files)	#we know how many spectra to plot
	#create lists of default options:
	colors=(nbspec*['blue','c','0.6','black'])[:nbspec]
	rvs=[0 for f in files]
	labels=['' for i,f in enumerate(files)]
	exts=[0 for f in files]
	columns=[1 for f in files]
	global labelOn; labelOn=False
	#read the list 'options' to see replace some of these default options:
	for i,opt in enumerate(options):
		if ('color' in opt):
			theColor=opt.split('=')[1]
			theIndex=int(opt.split('=')[0].replace('-color',''))
			colors[theIndex]=theColor
		elif ('rv' in opt):
			theRv=opt.split('=')[1]
			theIndex=int(opt.split('=')[0].replace('-rv',''))
			rvs[theIndex]=float(theRv)
		elif ('label' in opt):
			labelOn=True
			theLabel=opt.split('=')[1]
			theIndex=int(opt.split('=')[0].replace('-label',''))
			labels[theIndex]=theLabel
		elif ('ext' in opt):
			theExt=opt.split('=')[1]
			theIndex=int(opt.split('=')[0].replace('-ext',''))
			labels[theIndex]=int(theExt)
		elif ('column' in opt):
			theColumn=opt.split('=')[1]
			theIndex=int(opt.split('=')[0].replace('-column',''))
			columns[theIndex]=int(theColumn)

	#read the files one by one:
	for i,fileName in enumerate(files):
		if fileName[-5:]=='.fits':
			if noPyfits==False:
				#Read the fits file:
				hdulist = pyfits.open(fileName)
				#hdulist.info()			#displays info about the content of the file
								#(what we use for Daospec has only ONE extension)
				#print hdulist[0].header	#to print the whole header!
				wave_base = hdulist[exts[i]].header['CRVAL1']	# Angstrom
				try:
					wave_step = hdulist[exts[i]].header['CD1_1']	# Angstrom
				except:
					wave_step = hdulist[exts[i]].header['CDELT1']	# Angstrom
				flux = hdulist[exts[i]].data
				waveobs = np.arange(wave_base, wave_base+len(flux)*wave_step, wave_step)
				if len(waveobs) == len(flux) + 1:
					waveobs = waveobs[:-1]
				waveobs = [el*(1+(rvs[i]/c)) for el in waveobs]
				hdulist.close()
				wave_step1=wave_step
			else:
				print fileName,'cannot be read because PyFITS was not found.'
		else:
			#read the wavelength and flux from input txt:
			linesToRemove=1 #removes the first line because it is usually a header
			arr=[]
			try:
				inp = open(fileName,"r")
				lines = inp.readlines()[linesToRemove:]
				for line in lines:
					numbers = line.split()
					arr.append(numbers)
				waveobs = [float(el[0]) for el in arr]
				wave_step1 = (waveobs[-1]-waveobs[0])/(len(waveobs)-1)
				waveobs = [el*(1+(rvs[i]/c)) for el in waveobs]
				flux = [float(el[columns[i]]) for el in arr]
			except:
				#In case nothing was read at all:
				flux=[]
				waveobs=[]
				labels[i]=''
				print 'Problem reading the spectrum '+fileName




		#########          Plot it:
		plt.plot(waveobs, flux, color=colors[i], label=labels[i])








############ MAIN BODY:
#if the argument starts with a "@" then we batch-plot:
if len(sys.argv)==2 and sys.argv[1][0]=='@':
	instructionFile=sys.argv[1][1:]
	content=[]
	try:
		inp = open(instructionFile,"r")
		lines = inp.readlines()
		for line in lines:
			argz = line.split()
			content.append(argz)
	except:
		print 'Problem reading file',instructionFile
	i=0
	while 1>0:
		argumentsList=[sys.argv[0]]+content[i]
		plt.ion()
		plt.clf()
		actualseefits(argumentsList)
		plt.draw()
		spam=raw_input('Press ENTER for next plot.')
		if spam=='b':
			i=i-1	#to step back
		elif spam=='q':
			sys.exit()
		else:
			i=i+1	#to step forward
		if i==len(content):
			raw_input('Finished plotting.')
			sys.exit()
		
#if it doesn't start with a "@" then just a normal plot:
else:
	plt.figure(1)
	actualseefits(sys.argv)
	plt.ylabel('flux')
	plt.xlabel('wavelength')
	if labelOn==True:
		plt.legend()
	plt.show()
