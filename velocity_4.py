from starter2 import *
import tools.shot as shot
reload(shot)
import tools.tabletool as tabletool
reload(tabletool)

#names = ['r60','r120', 'r0', 's90','s120']
names = ['r60','r120', 'r0']
#names=['r60']

#checking for existence of things helps reduce work
if 'clobber' not in dir():
    clobber=False

#if a shot is in "device", it will not get remade.
if 'devices' not in dir() or clobber:
    devices={}

#lightmodel converts the image into usable density.
#Details are found in physical_values.py
lightmodel = {'r0':0,'r60':0,'r120':0,'s120':3,'s90':3}
#number of pixels for smoothing
smooth=3

#The following y-lines are used
y_cut = {'r0':[200,400], 'r60':[200,400], 'r120':[200,400], 's90':[450,650], 's120':[200,400]}

#A fine alignment is performed based on the voids in the preshock region.
#Both position and density are adjusted.  This region, in pixesls, is used to align.
bump_range={'r0':[175,300], 'r60':[175,300], 'r120':[175,300]}

#Region used for the shock front
vel_cut = {'r120':[325,575], 'r60':[325,550], 'r0':[410,600]}

#Points of the shock surround the foot and the peak of the shock.  For the shock jump.
shock_points = {'r60':[325,375, 650], 'r0':[380,500,700], 'r120':[325,375,650]}

#The post shock region position
post_shock_region = {'r0':[550,750], 'r60':[475,675], 'r120':[475,675]}

for name in names:
    if name not in devices:
        tmp=shot.device(name, lines=y_cut[name], lightmodel=lightmodel[name],smooth=smooth)
        devices[name]=tmp
        #tmp.image_density1(fname = 'image_shot_%s'%name)
        
        if name in ['s90','s120']:
            print("Analysis on %s incomplete, need fine adjust (bumper)"%name)
            continue

        x_off=None
        if name == 'r0':
            shift_60 = devices['r60'].shift_x
            shift_120 = devices['r120'].shift_x
            x_off = 0.5*(shift_60+shift_120)
        #Supply a file name to trigger a plot.
        tmp.bumper(bump_range[name],fix_shift_x=x_off    )#     ,fname='bumper_%s.pdf'%name)
        tmp.get_velocity( vel_cut[name]                  )#     ,fname = 'velocity_%s.pdf'%name)
        #tmp.get_velocity2( vel_cut[name]                      ,fname = 'velocity2_%s.pdf'%name)
        continue
        tmp.get_csound(mean_density=shock_points[name]        , fname='csound_%s.pdf'%name)
        tmp.get_sigma_rho(post_shock_region[name]             , fname = 'sigma_rho_%s.pdf'%name)
        #tmp.get_atwood()
        #tmp.get_sigma_v( )
        #print("Atwood %0.2f sigma_v %0.2f Ms %0.2f"%(tmp.atwood_number, tmp.sigma_v, tmp.sigma_v/tmp.cs))
        #print("mean rho %0.2f sigma_b %0.2f"%(tmp.mean_rho, tmp.sigma_B))

if 0:
    tabletool.table(devices, fname='table2.tex')
if 0:
    fig,ax=plt.subplots(1,1)
    for name in devices:
        dev = devices[name]
        ax.scatter( dev.sigma_v/dev.cs, dev.sigma_B, label=name)
    ax.legend(loc=0)
    ax.set(xlabel=r'$\sigma_v/c_s$', ylabel=r'$\sigma_B$')
    fig.savefig('%s/sigma_v_sigma_rho.pdf'%plot_dir)
