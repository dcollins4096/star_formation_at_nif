

from starter2 import *
import regions.regions_coarse as regions_coarse
import regions.regions_align as regions_align
reload(regions_coarse)

trimmed = regions_coarse.do_coarse_trim()
