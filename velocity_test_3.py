
from starter2 import *

import tools.enzo_gitter as eg
reload(eg)
plt.close('all')

reverse=False
particles=False
if 0:
    template = "simulations/advection/DD%04d/data%04d"
    name = 'advection'
    frame0=2
    frame1=3
    particles=True

if 1:
    template = "simulations/snowplow/DD%04d/data%04d"
    name = 'snowplow'
    frame0=259
    frame1=300
    reverse=True

if 0:
    name = 'bw'
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
    p_0 = gitter(frame0,'pressure')[sl]
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
    p_1 = gitter(frame1,'pressure')[sl]
    x_1   = gitter(frame1,'x')
    t_1   = gitter.ds['InitialTime']
    if particles:
        ad1   = gitter.ds.all_data()
        in_1  = ad1['particle_index']
        order = np.argsort(in_1)
        px_1  = ad1['particle_position_x'][order]
    delta_t = t_1-t_0
    delta_x = x_0[1]-x_0[0]
    delta_rho = rho_1-rho_0
    cs_y = gitter(frame0,'sound_speed')
    gamma = gitter.ds['Gamma']
    cs_0 = np.sqrt(gamma*p_0/rho_0)
    cs_1 = np.sqrt(gamma*p_1/rho_1)

import tools.velocity as vel
reload(vel)
fig,axes=plt.subplots(1,3, figsize=(12,4))
ax0=axes[0];ax1=axes[1]; ax2=axes[2]


sl=slice(None)
x_win,vvv1 = vel.vel1(rho_1[sl],rho_0[sl],dt=delta_t,dx=delta_x, ax=ax0)
x_won,vvv2 = vel.vel2(rho_1[sl],rho_0[sl],dt=delta_t,dx=delta_x)#, ax=ax0)
ax0.plot(x_0,rho_0,c='r')
ax0.plot(x_1,rho_1,c='g')
axt = ax0.twinx()
ok = np.abs(delta_rho)>1e-2 #hard coded for the snowplow
#ok = slice(None)
axt.plot( x_win[ok], vvv1[ok], c='purple')
axt.plot( x_won[ok], vvv2.v[ok], c='orange')
#axt.plot( x_0, 0.5*(vel_1+vel_0), c='pink')

if 1:
    shock_start = 0.75
    shock_end = 0.82
    e1=np.where(x_0<=shock_start)[0].max()
    e2=np.where(x_0<=shock_end)[0].max()
    shock_sl = slice(e1,e2)
    ax0.axvline(shock_start)
    ax0.axvline(shock_end)

    vs = vvv1[shock_sl]
    vg = vvv2[shock_sl]
    c_guess = np.sqrt( vs**2 - 0.5*vs*vg.v*(gamma+1))
    c_act = cs_0[shock_sl]
    mm = min([c_guess.min(),c_act.v.min()])
    xx = max([c_guess.max(),c_act.v.max()])
    #ax1.scatter(vs,vg)
    ax1.scatter(x_0[shock_sl],c_guess/c_act)
    ax1.plot(x_0[shock_sl],rho_0[shock_sl])
    ax1.plot(x_0[shock_sl],rho_1[shock_sl])
    #ax1.scatter(c_act, c_guess)
    #ax1.plot([mm,xx],[mm,xx])

if 0:
    import tools.equal_probability_binner as epb
    bins = np.linspace( vvv1[ok2].min(), vvv1[ok2].max(),16)
    ax1.hist(vvv1[ok],bins=bins, alpha=0.5)
    bins = np.linspace( -vvv2[ok2].max(), -vvv2[ok2].min(),16)
    ax1.hist(-vvv2[ok],bins=bins, alpha=0.5)

#pdf2,cen2,wid2=epb.equal_prob(np.log(vvv1), 16, ax=ax1)
#ax1.set(yscale='linear')


if 0:
    drho_dt = ((rho_1-rho_0)/delta_t)[1:-1]
    rhom = 0.5*(rho_1+rho_0)
    xcen = x_0[1:-1]
    drho_dx = (rhom[2:]-rhom[:-2])/(2*delta_x)
    c = np.zeros_like(drho_dt)
    ok = np.abs(drho_dx)>1e-10
    c[ok] = -drho_dt[ok]/drho_dx[ok]
    axt.plot( xcen[ok], c[ok], c='c')
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
fig.savefig('%s/vtest2_%s'%(plot_dir, name))
