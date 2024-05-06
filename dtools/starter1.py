#Main python modules.
#All external packages, none of my stuff.
import matplotlib
matplotlib.use('Agg')  #necessary some times.
import matplotlib.pyplot as plt
import matplotlib as mpl
import pylab
import h5py
import numpy as np
nar = np.array
import scipy
import astropy.io.fits as pyfits
import unyt

#need the reload module
import platform
ver=platform.python_version()
python_version=  int(ver[0])
if python_version == 3:
    from importlib import reload


#system things I want.
import sys
import re
import copy
import os
import copy
import warnings
import math
import pdb
import time
import glob
from collections import defaultdict

