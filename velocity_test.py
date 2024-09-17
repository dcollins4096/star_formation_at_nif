
from starter2 import *

import tools.enzo_gitter as eg
reload(eg)
plt.close('all')

reverse=False

if 0:
    template = "simulations/advection/DD%04d/data%04d"
    frame0=2
    frame1=3
if 1:
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


if 1:
    dx = x_0[1:]-x_0[:-1]
    dt = gitter.ds.arr(t_1-t_0,'s')
    drho_dt = (rho_1 - rho_0)/dt
    m = 0.5*(rho_1+rho_0)
    also_v = np.cumsum(-drho_dt*dx[0])/m 
    also_x = x_0
    vel_h = 0.5*(vel_0+vel_1)




fig,axes=plt.subplots(1,3,figsize=(8,4))
#ax0=axes[0][0];ax1=axes[0][1] 
#ax2=axes[1][0]#;ax3=axes[1][1]
ax0=axes[0];ax1=axes[1]; ax2=axes[2]
ax0.plot(x_0,rho_0.v)
ax0.plot(x_1,rho_1.v)
ax0.set(ylabel='rho')

ax1.plot(x_0,vel_0.v, label='v(t=0)',c='r')
#ax1.plot(x_1,vel_1.v, label='v(t=1)',c='g')

#ax1.plot(x_0,vel_h, label='v(t=1/2)',c='b')
ax1.plot(also_x,also_v,label='new',c='k')
ax1.legend(loc=0)
ax1.set(ylabel='vx')

#wtf
momentum = (rho_0*vel_0)

dpdx = (momentum[2:]-momentum[:-2])/(2*dx[0])
#ax2.plot(x_0[1:-1], -dpdx,c='r')
#ax2.plot(x_0, drho_dt, c='g')
rv = np.cumsum(dpdx*dx[0])
rv2 = np.cumsum(drho_dt*dx[0])
#ax2.plot(x_0[1:-1], rv, c='r')
ax2.plot(x_0, momentum, c='b')
ax2.plot(x_0, -rv2, c='g')
#axtwin=ax2.twinx()
ax2.plot(x_0,momentum+rv2)
#ax2.plot(x_0, momentum)
#ax2.plot(x_0,momentum/rho_0, c='r')
#ax2.plot(x_0,-rv2/rho_0,c='g')

fig.tight_layout()
fig.savefig('%s/velocity_from_density'%plot_dir)
