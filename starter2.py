
#starter 1 imports packages not specific to this project.

#Please install the following packages: pip works great.
#matplotlib
#numpy
#scipy
#astropy
#h5py
#unyt
#tifffile

from dtools.starter1 import *
import dtools.davetools as dt

from tifffile import imread
import tools.equal_probability_binner as epb
import regions.tiff_poker as TP

plot_dir = os.environ['HOME']+"/plots"
#plot_dir = "./plots_to_sort"
