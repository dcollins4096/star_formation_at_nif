

#
# Align the fiducial marks between pairs of images
#

from starter2 import *

import regions.aligner as aligner
reload(aligner)

def align_regions(trimmed, align_fname):

    TA_120_1, TA_120_2=aligner.align(trimmed['r120_t1'], trimmed['r120_t2'], name_base='r120', xrange=[290,350])
    TA_60_1, TA_60_2=aligner.align(trimmed['r60_t1'], trimmed['r60_t2'], name_base='r60', xrange=[300,340], x_override=0)
    TA_0_1, TA_0_2=aligner.align(trimmed['r0_t1'], trimmed['r0_t2'], name_base='r0', xrange=[300,330])

    fptr=h5py.File(align_fname,'w')
    fptr['r120_t1'] = TA_120_1
    fptr['r120_t2'] = TA_120_2
    fptr['r60_t1'] = TA_60_1
    fptr['r60_t2'] = TA_60_2
    fptr['r0_t1'] = TA_0_1
    fptr['r0_t2'] = TA_0_2
    fptr.close()


