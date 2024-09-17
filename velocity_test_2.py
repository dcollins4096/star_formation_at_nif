
from starter2 import *

import tools.enzo_gitter as eg
reload(eg)
plt.close('all')

reverse=False

if 1:
    template = "simulations/advection/DD%04d/data%04d"
    frame0=2
    frame1=3
if 0:
    template = "simulations/snowplow/DD%04d/data%04d"
    frame0=299
    frame1=300
    reverse=True

if 0:
    template = "simulations/BW/DD%04d/data%04d"
    frame0=30
    frame1=31
    reverse=False

if 'rho_0' not in dir():
    gitter = eg.gitter2(template)
    sl = slice(None)
    velrev=1
    if reverse:
        sl = slice(-1,None,-1)
        velrev=-1
    rho_0 = gitter(frame0,'density')[sl]
    x_0   = gitter(frame0,'x')
    vel_0 = velrev*gitter(frame0,'x-velocity')[sl]
    t_0   = gitter.ds['InitialTime']
    rho_1 = gitter(frame1,'density')[sl]
    vel_1 = velrev*gitter(frame1,'x-velocity')[sl]
    x_1   = gitter(frame1,'x')
    t_1   = gitter.ds['InitialTime']
    delta_t = t_1-t_0
    delta_x = x_0[1]-x_0[0]


import tools.velocity as vel
reload(vel)
fig,axes=plt.subplots(1,2)
ax0=axes[0];ax1=axes[1]

vvv1 = vel.vel1(rho_1,rho_0,dt=delta_t,dx=delta_x, ax=ax0)
import tools.equal_probability_binner as epb
pdf2,cen2,wid2=epb.equal_prob(vvv1, 16, ax=ax1)
fig.tight_layout()
fig.savefig('%s/vtest2'%plot_dir)
