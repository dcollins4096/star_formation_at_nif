from starter2 import *
import tools.shot as shot
reload(shot)

names = ['r60','r120', 'r0']
#names=['r60']

if 'devices' not in dir() :
    devices={}

#The following y-lines are used
y_cut = {'r0':[200,400], 'r60':[200,400], 'r120':[200,400]}

#A fine alignment is performed based on the voids in the preshock region.
#Both position and density are adjusted.  This region, in pixesls, is used to align.
bump_range={'r0':[175,300], 'r60':[175,300], 'r120':[175,300]}

#Region used for the shock front
vel_cut = {'r120':[325,575], 'r60':[325,550], 'r0':[410,600]}

#Points of the shock surround the foot and the peak of the shock.  For the shock jump.
shock_points = {'r60':[325,375, 650], 'r0':[380,500,700], 'r120':[325,375,650]}

#The post shock region position
post_shock_region = {'r0':[550,750], 'r60':[475,675], 'r120':[475,675]}

#lightmodel converts the image into usable density.
#Details are found in physical_values.py
lightmodel=0
#number of pixels for smoothing
smooth=3
for name in names:
    if name not in devices:
        tmp=shot.device(name, lines=y_cut[name], lightmodel=lightmodel,smooth=smooth)
        tmp.image_density(fname = 'image_shot_%s'%name)
        x_off=None
        if name == 'r0':
            shift_60 = devices['r60'].shift_x
            shift_120 = devices['r120'].shift_x
            x_off = 0.5*(shift_60+shift_120)
        tmp.bumper(bump_range[name],fix_shift_x=x_off,fname='bumper_%s.pdf'%name)
        tmp.compute_velocity( vel_cut[name],fname = 'velocity_%s'%name)
        tmp.sigma_rho(post_shock_region[name], fname = 'sigma_rho_%s'%name)
        tmp.csound(mean_density=shock_points[name], fname='csound_%s'%name)
        tmp.atwood()
        tmp.sigma_v( )
        devices[name]=tmp

if 1:
    fig,ax=plt.subplots(1,1)
    for name in devices:
        dev = devices[name]
        ax.scatter( dev.sigma_v/dev.cs, dev.sigma_B, label=name)
    ax.legend(loc=0)
    fig.savefig('%s/sigma_v_sigma_rho.pdf'%plot_dir)
