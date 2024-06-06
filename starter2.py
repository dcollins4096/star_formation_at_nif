

#Please install the following packages: pip works great.
#matplotlib
#numpy
#scipy
#astropy
#h5py
#unyt
#tifffile

#Load packages not specific to this project.
import os
import sys
if not os.path.exists('dtools/starter1.py'):
    print("Need to get the submodule")
    print("git submodule init")
    print("git submodule update")
    sys.exit(-1)
from dtools.starter1 import *
import dtools.davetools as dt
from scipy.ndimage import gaussian_filter

#make this a directory that makes sense for you
plot_dir = os.environ['HOME']+"/plots"
#plot_dir = "./plots_to_sort"

#project specific imports
import tools.equal_probability_binner as epb
import tools.power_spectrum as ps
import tools.horizontal_distance as horz
