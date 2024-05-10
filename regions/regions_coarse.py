#
# Rough cut the data from the origin TIFF.
#

from starter2 import *
import regions.tiff_poker as TP
reload(TP)

def raw_image(which, fname):
    #plot the regions
    V = TP.viewer(which=which)
    V.guess_scale()
    #V.image(vmin=V.minmax[0],vmax=V.minmax[1])
    V.image1('%s/%s'%(plot_dir,fname))



def do_coarse_trim():
#not the most elegant of code.
#a,b,c,d are left, right, bottom, and top in pixels.
    trim = {}
    first_round = [0,1,2,3,4,5]
    second_round = [8,9,10,11]
    both = first_round + second_round
    for DO in both:
        if DO==0:
            which=0
            a=550;b=1800;c=300;d=1000
            na=0;nb=600;nc=0;nd=1000
            section = 'r0_t1'
            zero = 'Negative'

        if DO==1:
            which=0
            a=100;b=1700;c=1100;d=1900
            na=0;nb=600;nc=0;nd=1000
            section = 'r0_t2'
            zero = 'Negative'

        if DO==2:
            which=1
            section = 'r60_t1'
            a=550;b=1800;c=300;d=1000
            na=50;nb=600;nc=0;nd=1000
            zero = 'Negative'

        if DO==3:
            which=1
            section = 'r60_t2'
            a=100;b=1700;c=1100;d=1900
            na=50;nb=600;nc=0;nd=1000
            zero = 'Negative'


        if DO==4:
            which=2
            section = 'r120_t1'
            a=550;b=1800;c=300;d=1100
            na=50;nb=600;nc=0;nd=1000
            zero = 'Negative'

        if DO==5:
            which=2
            section = 'r120_t2'
            a=100;b=1700;c=1100;d=2200
            na=50;nb=600;nc=0;nd=1000
            zero = 'Negative'

        #day 2
        if DO==6:
            #has the huge arc
            which=3
            section = 's0_t1'
            a=375;b=1200;c=50;d=950
            na=600;nb=1200;nc=950;nd=1400
            zero = 'Positive'

        if DO==7:
            #the other frame from the arc shot.  
            #values not right.
            which=3
            section = 's0_t2'
            a=375;b=1200;c=50;d=950
            na=600;nb=1200;nc=950;nd=1400
            zero = 'Positive'

        if DO==8:
            which=4
            section = 's90_t1'
            a=275;b=1200;c=60;d=950
            na=600;nb=1200;nc=950;nd=1400
            zero = 'Positive'
        if DO==9:
            which=4
            section = 's90_t2'
            a=1200;b=2100;c=100;d=1700
            na=600;nb=1200;nc=950;nd=1400
            zero = 'Positive'
        if DO==10:
            which=5
            section = 's120_t1'
            a=275;b=1200;c=60;d=950
            na=600;nb=1200;nc=950;nd=1400
            zero = 'Positive'
        if DO==11:
            which=5
            section = 's120_t2'
            a=1200;b=2100;c=100;d=1700
            na=600;nb=1200;nc=950;nd=1400
            zero = 'Positive'

        print("extract ", section)
        if 1:
            #plot the regions
            V = TP.viewer(which=which)
            V.guess_scale(section=section)
        if 1:
            #V.image(vmin=V.minmax[0],vmax=V.minmax[1])
            V.image1('%s/%s_full_image'%(plot_dir,section))

        if 1:
            #Get the noise level for trimming
            x,y,noise=V.xtract(a=na,b=nb,c=nc,d=nd)
            vmax = noise.max()
            sigma_n = noise.std()
            x,y,noise=V.xtract_and_image(a=na,b=nb,c=nc,d=nd,vmin=0,vmax=vmax, fname='%s/%s_noises'%(plot_dir,section),zero=zero)

        if 0:
            #The brigt one is hard to deal with.
            x,y,rough=V.xtract_and_image(a=a,b=b,c=c,d=d,vmin=V.minmax[0],vmax=V.minmax[1], fname='%s/%s_rough'%(plot_dir,section),zero=False)
            trim[section] = TP.trimmer2(rough, fname="%s/%s_trim"%(plot_dir,section),zero=zero)
        if 1:
            #auto-trim based on the noise level.
            x,y,rough=V.xtract_and_image(a=a,b=b,c=c,d=d,vmin=V.minmax[0],vmax=V.minmax[1], fname='%s/%s_rough'%(plot_dir,section),zero=False)
            trim[section] = TP.trimmer(rough, sigma_n=2*sigma_n,fname='%s/%s_trim'%(plot_dir,section), vmin=V.minmax[0],vmax=V.minmax[1], zero=zero)
            if section.startswith('s'):
                trim[section]=trim[section].transpose()[:,::-1]
    return trim

def save_trim(trim,fname):
    fptr=h5py.File(fname,'w')
    for section in trim:
        fptr[section]=trim[section]
    fptr.close()
    print('saved file',fname)

def read(fname):
    fptr=h5py.File(fname,'r')
    trim={}
    for section in fptr:
        trim[section]=fptr[section][()]
    fptr.close()
    return trim

