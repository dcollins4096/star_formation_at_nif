
from starter2 import *

import tools.enzo_gitter as eg
reload(eg)
plt.close('all')

reverse=False
particles=False
if 1:
    template = "simulations/advection/DD%04d/data%04d"
    frame0=2
    frame1=3
    particles=True
if 0:
    template = "simulations/snowplow/DD%04d/data%04d"
    frame0=259
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
    if particles:
        ad0   = gitter.ds.all_data()
        in_0  = ad0['particle_index']
        order = np.argsort(in_0)
        px_0  = ad0['particle_position_x'][order]
    rho_1 = gitter(frame1,'density')[sl]
    vel_1 = velrev*gitter(frame1,'x-velocity')[sl]
    x_1   = gitter(frame1,'x')
    t_1   = gitter.ds['InitialTime']
    if particles:
        ad1   = gitter.ds.all_data()
        in_1  = ad1['particle_index']
        order = np.argsort(in_1)
        px_1  = ad1['particle_position_x'][order]
    delta_t = t_1-t_0
    delta_x = x_0[1]-x_0[0]

if 0:
    rho_0 = 0*x_0
    rho_1 = 0*x_0
    n1=30
    n2=20
    rho_0[n1:] =rho_0[n1].v+ x_0[n1:].v-x_0[n1].v
    rho_1[n2:] =rho_0[n2].v+ x_0[n2:].v-x_0[n2].v


import tools.velocity as vel
reload(vel)
fig,axes=plt.subplots(1,3, figsize=(12,4))
ax0=axes[0];ax1=axes[1]; ax2=axes[2]


sl = slice(0,40)
sl = slice(None)
x_win,vvv1 = vel.vel1(rho_1[sl],rho_0[sl],dt=delta_t,dx=delta_x)#, ax=ax0)
x_won,vvv2 = vel.vel2(rho_1[sl],rho_0[sl],dt=delta_t,dx=delta_x)#, ax=ax0)
ax0.plot(x_0,rho_0,c='r')
ax0.plot(x_1,rho_1,c='g')
axt = ax0.twinx()
axt.plot( x_win, -vvv1, c='purple')
axt.plot( x_won, -vvv2.v+vel_0[0].v, c='orange')
axt.plot( x_0, 0.5*(vel_1+vel_0), c='pink')
import tools.equal_probability_binner as epb
bins = np.linspace( vvv1.min(), vvv1.max(),16)
ax1.hist(vvv1,bins=bins)
#pdf2,cen2,wid2=epb.equal_prob(np.log(vvv1), 16, ax=ax1)
#ax1.set(yscale='linear')


if 1:
    drho_dt = ((rho_1-rho_0)/delta_t)[1:-1]
    rhom = 0.5*(rho_1+rho_0)
    xcen = x_0[1:-1]
    drho_dx = (rhom[2:]-rhom[:-2])/(2*delta_x)
    c = np.zeros_like(drho_dt)
    ok = np.abs(drho_dx)>1e-10
    c[ok] = -drho_dt[ok]/drho_dx[ok]
    axt.plot( xcen[ok], c[ok], c='c')
    cs = gitter(frame0,'sound_speed')
    ok2 = np.abs(vvv2) > 1e-12
    c2 = np.zeros_like(vvv2)
    c2[ok2] = (cs.v**2+vvv2.v**2)[ok2]/vvv2[ok2]
    ax2.plot( np.abs(c))
    ax2.plot( np.abs(vvv1))
    #ax2.plot( np.abs(c2))
    #ax2.set(yscale='log')


if particles and False:
    ax0.scatter(px_0, [1]*px_0.size, c='r')
    ax1.scatter(px_1, [1]*px_1.size, c='g')

fig.tight_layout()
fig.savefig('%s/vtest2'%plot_dir)
