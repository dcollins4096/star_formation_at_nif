
from starter2 import *

import tools.enzo_gitter as eg
reload(eg)
plt.close('all')

template = "simulations/advection/DD%04d/data%04d"
frame0=2
frame1=3
if 'rho_0' not in dir():
    gitter = eg.gitter(template)
    rho_0 = gitter(frame0,'density')
    x_0   = gitter(frame0,'x')
    vel_0 = gitter(frame0,'x-velocity')
    t_0   = gitter.ds['InitialTime']
    rho_1 = gitter(frame1,'density')
    vel_1 = gitter(frame1,'x-velocity')
    x_1   = gitter(frame1,'x')
    t_1   = gitter.ds['InitialTime']


dx = x_0[1:]-x_0[:-1]
dt = gitter.ds.arr(t_1-t_0,'s')
drho_dt = (rho_1 - rho_0)/dt
m = rho_1
also_v = np.cumsum(-drho_dt*dx[0])/m + vel_0[0]
also_x = x_0
vel_h = 0.5*(vel_0+vel_1)



fig,axes=plt.subplots(1,2,figsize=(8,4))
#ax0=axes[0][0];ax1=axes[0][1] 
#ax2=axes[1][0]#;ax3=axes[1][1]
ax0=axes[0];ax1=axes[1]#; ax2=axes[2]
ax0.plot(x_0,rho_0.v-1)
ax0.plot(x_1,rho_1.v-1)
ax0.set(ylabel='rho')

ax1.plot(x_0,vel_0.v, label='v(t=0)')
ax1.plot(x_1,vel_1.v, label='v(t=1)')
ax1.set(ylabel='vx')

ax1.plot(x_0,vel_h, label='v(t=1/2)')
ax1.plot(also_x,also_v,label='RHO',c='k')

fig.tight_layout()
fig.savefig('%s/velocity_from_density'%plot_dir)
