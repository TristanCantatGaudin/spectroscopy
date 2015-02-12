#!/usr/bin/env python
#Tristan Cantat-Gaudin, 28/10/2014.
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
#
#		-column0=1		for ascii spectra, which column contains the flux
#					(default is 1, and the wavelength has to be column 0)
#
#		-dark			for a black/red/orange/yellow color scheme
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
#			Type "!comment" to store comments in a log file
#			Type b to go back to the previous line.
#			Type q to exit.
#
#	On the plot window: keep "v" key pressed and click on H-alpha or H-beta
#	line to estimate the radial velocity. (non-relativistic Doppler)
#
#			    keep "f" pressed and click points to follow a line profile
#			    then keep "d" pressed and click anywhere to fit a gaussian profile
#			    keep "x" pressed and click anywhere to erase all the fitted profiles
#			    (this only works when plotting FITS spectra!)
#
#	EXAMPLES:
#		./seefits.py file.fits
#
#		./seefits.py file.fits[4]
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


dark = '-dark' in sys.argv
if '-steps' in sys.argv:
	drawstyle='steps-mid'
else:
	drawstyle=''



#--------------------------- FUNCTIONS ---------------------------------------
#Defines what happens when a part of the plot is clicked on:
global clicksForGaussian; clicksForGaussian=[]
global plotElementsToDelete; plotElementsToDelete=[]
global wave_step
####################
def on_click(event):
	global clicksForGaussian, plotElementsToDelete, wave_step
	if event.key=='v':
		wclick = event.xdata
		line,vel = wav2radvel(wclick)
		print vel,'km/s (using',str(line)+')'

	if event.key=='f': #keep f key pressed and click to follow a line profile
		clicksForGaussian.append( [event.xdata,event.ydata] )
		elementToDelete, = ax.plot(event.xdata,event.ydata,'ro')
		plotElementsToDelete.append( elementToDelete )
		#print clicksForGaussian
	if event.key=='d': #click somewhere while pressing d to fit a gaussian to the defined set of points
		if len(clicksForGaussian)==0:
			print 'Keep f pressed and click points, then keep d pressed and click anywhere to fit a gaussian.'
		else:
			from scipy.optimize import curve_fit
			xtofit,ytofit=zip(*clicksForGaussian)
			def gaus(x,a,x0,sigma,y0):
				return y0-a*np.exp(-(x-x0)**2/(2*sigma**2))
			popt,pcov = curve_fit(gaus,xtofit,ytofit,p0=[max(ytofit)-min(ytofit),np.mean(xtofit),np.std(xtofit),max(ytofit)])
			xtoplot=np.linspace(min(xtofit),max(xtofit),100)
			elementToDelete, = ax.plot(xtoplot,gaus(xtoplot,*popt),'r-',lw=2)
			plotElementsToDelete.append( elementToDelete )
			EW=popt[0]*popt[2]*np.sqrt(2*np.pi)/max(ytofit) #a*sigma*sqrt(2pi), normalising to continuum level!
			print 'WL=%8.3f  FWHM=%5.3fmA =%3.1fpx  EW=%5.3fmA' % (popt[1],10*2.3548*popt[2],2.35482*popt[2]/wave_step,1000*EW)
			clicksForGaussian=[]
	if event.key=='x': #click somwhere while pressing x to delete this line
		for el in plotElementsToDelete:
			ax.lines.remove(el)
		plotElementsToDelete=[]
		print 'Cleared.'




def wav2radvel(wavobs,wavrest=0):
	'''
	Input: observed wavelength, rest wavelength
	Output: radial velocity
	If no rest wavelength given, it assumes you mean H-alpha or H-beta (whichever is closest).
	'''
	wavHa=6562.797
	wavHb=4861.323
	if wavrest==0:
		dist = wavobs - ((wavHa+wavHb)/2)
		if dist>=0:
			line = 'H-alpha'
			wavrest = wavHa
		else:
			line = 'H-beta'
			wavrest = wavHb
	radvel=299792.458*( (1.*wavobs/wavrest) - 1)
	return line,radvel







#function that takes the arguments, reads the files and plots:
def actualseefits(argumentsList):
	global wave_step
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
	if dark==True:
		colors=(nbspec*['white','r','#FF7F00','#FFFF00'])[:nbspec]
	else:
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
		elif ('column' in opt):
			theColumn=opt.split('=')[1]
			theIndex=int(opt.split('=')[0].replace('-column',''))
			columns[theIndex]=int(theColumn)

	#obtain list of extensions to read:
	for i,fileName in enumerate(files):
		#find in requested a specific extension, ext=0 if not:
		try:	
			toto=fileName[-3:]
			if toto[0]=='[' and toto[-1]==']':
				exts[i]=int(toto[1])
				files[i]=fileName[:-3]
			else:
				exts[i]=0
		except:
			exts[i]=0

	#read the files one by one:
	for i,fileName in enumerate(files):
		if fileName[-5:]=='.fits':
			if noPyfits==False:
				#Read the fits file:
				hdulist = pyfits.open(fileName)
				#hdulist.info()			#displays info about the content of the file
								#(what we use for Daospec has only ONE extension)
				#print hdulist[0].header	#to print the whole header!

				try:
					pix_start = hdulist[exts[i]].header['CRPIX1']   # starting pixel
				except:
					pix_start = 1
				try:
					wave_step = hdulist[exts[i]].header['CD1_1']	# Angstrom
				except:
					wave_step = hdulist[exts[i]].header['CDELT1']	# Angstrom
				wave_base = hdulist[exts[i]].header['CRVAL1']	# Angstrom
				wave_base=wave_base+(1-pix_start)*wave_step
				#(necessary in case CRPIX1 is not 1)
				flux = hdulist[exts[i]].data
				waveobs = np.arange(wave_base, wave_base+len(flux)*wave_step, wave_step)
				#these ligns can solve issues with NumPy "arange":
				if len(waveobs) == len(flux) + 1:
					waveobs = waveobs[:-1]
				#
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
		plt.plot(waveobs, flux, color=colors[i], label=labels[i],drawstyle=drawstyle)
	return files[0]	#return the name of the first file







############ MAIN BODY:
#if the argument starts with a "@" then we batch-plot:
if len(sys.argv)==2 and sys.argv[1][0]=='@':
	logFileName=sys.argv[1][1:]+'.log'	# <-----! file in which to store the comments
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
		title = actualseefits(argumentsList)
		plt.draw()
		if i==len(content)-1:
			spam=raw_input(title+'\t  This is the last spectrum. ')
		else:
			spam=raw_input(title+'\t  Press ENTER for next plot. ')
		if spam=='b':
			i=i-1	#to step back
		elif spam=='q':
			sys.exit()
		elif ('!' in spam):
			if spam[0]=='!':	#to write comments
				logFile = open(logFileName, "a")
				comment=spam[1:]
				logFile.write(title+'\t'+comment+'\n')
				logFile.close()
		else:
			i=i+1	#to step forward

		if i==len(content):
			i=i-1
		
#if the argument doesn't start with a "@" then just a normal plot:
else:
	if dark==True:
		fig=plt.figure(1,facecolor='black', edgecolor='black')
		import matplotlib as mpl
		mpl.rc('axes',facecolor='black',edgecolor='white',labelcolor='white')
		mpl.rc('text',color='white')
		mpl.rc('xtick',color='white')
		mpl.rc('ytick',color='white')
		mpl.rc('savefig',facecolor='black',edgecolor='black')
	else:
		fig=plt.figure(1)
	title = actualseefits(sys.argv)
	ax=fig.add_subplot(111)
	plt.ylabel('flux')
	plt.xlabel('wavelength')
	plt.title(title)
	if labelOn==True:
		plt.legend()
	fig.canvas.mpl_connect('button_press_event', on_click) #to handle clicks!
	plt.show()
