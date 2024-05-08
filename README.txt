Prereqs:
--------
#Please install the following packages: pip works great.
#
matplotlib
numpy
scipy
astropy
h5py
unyt
tifffile


Directories/modules:
-------------------
data:    Put the raw tiff files in here.  Can I check those in?
dtools:  some cross platform tools
tools:   project specific tools
regions: contains tools to extract regions of the data.
plots:   a place for your plots.

starter2.py: loads things
trim_and_align.py: extract frames from raw tiff.
image_subregions: take image of preshock and noise regions for calibrating the image.
plot_shock.py:  plot the shock regions 
	        convert image to density
		small scale align (bumper)
		compute velocity from front motion
		compute sound speed from density jump
		compute density variance with Brunt method
		compute atwood number
		compute sigma_v
	       
