
from dtools.starter1 import *
import power_spectrum as ps

import scipy.stats
import equal_probability_binner as ep

from tifffile import imread
import tiff_poker as TP
reload(TP)
plt.close('all')

if 'TNA' not in dir():
    fname = 'p68_laser/TRIM_ALIGN.h5'
    fptr=h5py.File(fname,'r')
    TNA = {}
    for field in fptr:
        TNA[field]=fptr[field][()]
    fptr.close()

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
    for shot in zero_abcd:
        which = {'r0_t1':0,'r0_t2':0,'r60_t1':1,'r60_t2':1,'r120_t1':2,'r120_t2':2}[shot]
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
        which = {'r0_t1':0,'r0_t2':0,'r60_t1':1,'r60_t2':1,'r120_t1':2,'r120_t2':2}[shot]
        V = TP.viewer(which=which)
        a,b,c,d=zero_abcd[shot]
        fname = 'plots_to_sort/zero_%s'%shot
        X,Y,Z = V.xtract_and_image(a=a,b=b,c=c,d=d,vmin=None,vmax=None,fname=fname, zero=True)
        print(zero_region[shot]-Z)

def get_zero(shot):
    zero = zero_region[shot]
    hist, cen = ep.equal_prob(zero.flatten(), 16)
    zero_val = cen[np.argmax(hist)]
    return zero_val



def image_preshock():
    for shot in preshock_abcd:
        arr=TNA[shot]
        V = TP.viewer(arr=arr)
        a,b,c,d=preshock_abcd[shot]
        fname = 'plots_to_sort/preshock_%s'%shot
        X,Y,Z = V.xtract_and_image(a=a,b=b,c=c,d=d,vmin=None,vmax=None,fname=fname)


if 0:
    #denoise and make an image
    plot_dir='plots_to_sort'
    import get_foam as gf
    reload(gf)
    denoise_preshock={}

    for shot in preshock_region:
        NOI = gf.noisy(preshock_region[shot])
        NOI.plot_denoise(kmin=0,nmax=9,outname='%s/%s_denoise'%(plot_dir,shot))
        denoise_preshock[shot] = NOI.rhoback

                                         
