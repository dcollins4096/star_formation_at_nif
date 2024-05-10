

#
# Align the fiducial marks between pairs of images
#

from starter2 import *

import regions.aligner as aligner
reload(aligner)

def align_regions2(trimmed, align_fname=None, shot_list=None):
    shot_xrange={'r0':[300,330], 'r60':[300,340], 'r120':[290,350], 's120':[200,300], 's90':[200,400]}
    thresh = {'s120':100,'s90':100}
    fid_line = {'s120':50, 's90':50}
    shots={}

    preserve_x = {'r0':False,'r60':False,'r120':False,'s90':True,'s120':True}

    for name in shot_list:
        fname = "pre_align_%s"%name
        name1 = name+"_t1"
        name2 = name+"_t2"
        aligner.image_pre_align(trimmed[name1],trimmed[name2],fname)

    #TA_120_1, TA_120_2=aligner.align(trimmed['r120_t1'], trimmed['r120_t2'], name_base='r120', xrange=[290,350])
    #TA_60_1, TA_60_2=aligner.align(trimmed['r60_t1'], trimmed['r60_t2'], name_base='r60', xrange=[300,340], x_override=0)
    #TA_0_1, TA_0_2=aligner.align(trimmed['r0_t1'], trimmed['r0_t2'], name_base='r0', xrange=[300,330])
    for name in shot_list:
        shots[name] = aligner.align(trimmed[name+'_t1'], trimmed[name+'_t2'], name_base=name, xrange=shot_xrange[name], 
                                    thresh=thresh.get(name,250), fid_line=fid_line.get(name,600), preserve_x=preserve_x[name])

    fptr=h5py.File(align_fname,'w')
    for name in shot_list:
        fptr[name+'_t1']= shots[name][0]
        fptr[name+'_t2']= shots[name][1]
    fptr.close()
    print('save',align_fname)
def align_regions_old(trimmed, align_fname):

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



