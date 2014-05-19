Python Spectroscopy Scripts
===========================


Python is a very convenient language for scripting, plotting, manipulating data, reading and writing files. I share here scripts that I have written for my personal use (mainly for spectroscopy) and that may be useful to other people. Some of them can work with either ASCII or FITS data, so if you want to use them on FITS files you need to have PyFITS installed.


plotSpec.py
-----------
This can plot spectra contained in FITS files, or in ASCII tables. It is convenient because it can plot several spectra from different formats at the same time. More details about colors, radial velocity correction and reading format at the beginning of the file. Ex:

    ./plotSpec.py spectrum.fits spectrum_deg.fits -label0=original -label1=degraded

or:

    ./plotSpec.py arcturus.fits[2]


seehead.py
-----------
Shows the structure and header of a FITS file. Ex:

    ./seehead.py spectrum.fits

to show the list of data units (list of extensions) contained in a FITS file.

    ./seehead.py spectrum.fits 0

to print the header of the extension 0.


txt2fits.py
-----------
Convert a spectrum in an ASCII table to a fits file. NB: the FITS format requires a constant step in wavelength, so if your data is not sampled on a constant step there is an option to resample it (in its current version the script is quite slow at doing that, but works fine). More details about the options (especially the clean the spectrum, fill in gaps, avoid negative fluxes, cut spikes, add some radial velocity shift etc) in the file. Ex:

    ./txt2fits.py spectrum.txt rv=35 


fits2txt.py
-----------
Extract a FITS spectrum to a two-column ASCII table. Ex:

    ./fits2txt.py spectrum.fits

or:

    ./fits2txt.py spectrum.fits[3]


cutSpec.py
-----------
Get a slice of a FITS file spectrum, specifying a wavelength range. By default the extension 0 is read, but you can also specify which extension to read (and use this script to extract a single extension from a multi-ext. FITS). Ex:

    ./cutSpec.py largespectrum.fits 4800 5800 shorterspectrum.fits

or:

    ./cutSpec.py largespectrum[2].fits 4800 5800 shorterspectrum.fits


degrade.py
----------
Degrade the resolution of a spectrum. If your spectrum has a fixed resolution, it means the fw of the lines increases with the wavelength. Degrading the resolution is not just a convolution, because the width of the gaussian used for the convolution scales with the wavelength. This script reads an ASCII or FITS spectrum and writes the degraded one in a new file (of the same format as the input). Ex:

    ./degrade.py spectrum.fits 80000 47000

The degrade_f.py script is the same, but uses a Fortran subroutine for the inner loop, that is compiled with f2py. The subroutine is actually included in the script itself, and is written and compiled at the first execution.

