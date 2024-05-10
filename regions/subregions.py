
from starter2 import *

import regions.tiff_poker as TP
import regions.regions_coarse as regions_coarse

#get the Trim And Aligned values
if 'TNA' not in dir():
    trim_fname = "data/aligned.h5"
    TNA = regions_coarse.read(trim_fname)

if 'preshock_region' not in dir():
    preshock_region={}
    preshock_abcd={}
    #abcd = left, right, bottom, top
    preshock_abcd['r60_t1'] = [100,400,100,400]
    preshock_abcd['r60_t2'] = [100,400,100,400]

    preshock_abcd['r120_t1'] = [100,400,100,400]
    preshock_abcd['r120_t2'] = [100,400,100,400]

    preshock_abcd['r0_t1'] = [100,400,100,400]
    preshock_abcd['r0_t2'] = [100,400,100,400]

    preshock_abcd['s90_t1'] = [100,400,150,400]
    preshock_abcd['s90_t2'] = [100,400,150,400]
    preshock_abcd['s120_t1'] = [100,400,100,400]
    preshock_abcd['s120_t2'] = [100,400,100,400]
    for shot in preshock_abcd:
        arr=TNA[shot]
        V = TP.viewer(arr=arr)
        a,b,c,d=preshock_abcd[shot]
        X,Y,Z = V.xtract(a=a,b=b,c=c,d=d)
        preshock_region[shot] = Z

if 'zero_region' not in dir():
    zero_region={}
    zero_abcd={}
    #abcd = left, right, bottom, top
    zero_abcd['r60_t1'] = [100,600,100,750]
    zero_abcd['r60_t2'] =  [100,600,100,750]

    zero_abcd['r120_t1'] =  [100,600,100,750]
    zero_abcd['r120_t2'] =  [100,600,100,750]

    zero_abcd['r0_t1'] = [100,600,100,750]
    zero_abcd['r0_t2'] = [100,600,100,750]

    zero_abcd['s90_t1']  = [600,1200,950,1400]
    zero_abcd['s90_t2']  = [600,1200,950,1400]
    zero_abcd['s120_t1'] = [600,1200,950,1400]
    zero_abcd['s120_t2'] = [600,1200,950,1400]
    for shot in zero_abcd:
        which = {'r0_t1':0,'r0_t2':0,'r60_t1':1,'r60_t2':1,'r120_t1':2,'r120_t2':2,'s90_t1':4, 's90_t2':4, 's120_t1':5,'s120_t2':5}[shot]
        V = TP.viewer(which=which)
        a,b,c,d=zero_abcd[shot]
        fname = 'plots_to_sort/zero_%s'%shot
        X,Y,Z = V.xtract(a=a,b=b,c=c,d=d)
        #X,Y,Z = V.xtract_and_image(a=a,b=b,c=c,d=d,vmin=None,vmax=None,fname=fname, zero=True)
        zero_region[shot] = Z

def image_zero(shot_list = None):
    if shot_list is None:
        shot_list = list(zero_abcd.keys())
    for shot in shot_list:
        print('image zero ',shot)
        which = {'r0_t1':0,'r0_t2':0,'r60_t1':1,'r60_t2':1,'r120_t1':2,'r120_t2':2,'s90_t1':4, 's90_t2':4, 's120_t1':5,'s120_t2':5}[shot]
        V = TP.viewer(which=which)
        a,b,c,d=zero_abcd[shot]
        fname = '%s/zero_%s'%(plot_dir,shot)
        X,Y,Z = V.xtract_and_image(a=a,b=b,c=c,d=d,vmin=None,vmax=None,fname=fname, zero=True)

def get_zero(shot):
    zero = zero_region[shot]
    hist, cen = epb.equal_prob(zero.flatten(), 16)
    zero_val = cen[np.argmax(hist)]
    return zero_val

def image_preshock(shot_list=None):
    if shot_list is None:
        shot_list = list(preshock_abcd.keys())

    for shot in shot_list:
        print("Image",shot)
        arr=TNA[shot]
        V = TP.viewer(arr=arr)
        a,b,c,d=preshock_abcd[shot]
        fname = '%s/preshock_%s'%(plot_dir,shot)
        X,Y,Z = V.xtract_and_image(a=a,b=b,c=c,d=d,vmin=None,vmax=None,fname=fname, zero=False)


