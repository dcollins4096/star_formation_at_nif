from dtools.starter1 import *
import tools.horizontal_distance as horz
reload(horz)

def vel2(rho1, rho2, dx=1, dt=1, ax=None):
    x = np.arange(rho1.size)*dx
    drho_dt = (rho2-rho1)/dt
    m = 0.5*(rho2+rho1)
    v = np.cumsum(-drho_dt*dx)/m 
    return x, v

def vel1(rho1, rho2, dx=1, dt=1, ax=None):
    """get the velocity between two densities using lines at constant density"""
    I2 = horz.ho2(ya=rho1.v, yb=rho2.v)
    #I2 = horz.ho(rho1.v, rho2.v)
    dx2 = I2[:,1]-I2[:,0]
    velocity = dx2*dx/dt
    print('dx/dt', dx/dt)
    if 1:
        fig,axes=plt.subplots(1,1)
        ax0=axes
        #ax0.plot( rho1[nar(I2[:,0].astype('int'))])
        ax0.plot( (I2[:,1] - I2[:,0]))
        #ax0.plot( I2[:,0])
        fig.savefig('%s/vel_dumb'%plot_dir)
    if ax is not None:
        thex = np.arange(len(rho1))*dx
        #ax.plot(thex, rho2-rho1,c='r',marker='*')
        thex = np.arange(len(rho1))*dx
        ax.plot(thex, rho1,c='r')
        ax.plot(thex, rho2,c='g')
        for ii, x, y in I2:
            i=int(ii)
            xs = [i*dx, x*dx]
            ys = [rho1[i], y]
            ax.plot(xs, ys, "g--")
        ax2=ax.twinx()
        xxx = np.arange(velocity.size)*dx
        ax2.plot(xxx, velocity,c='purple')

    return I2[:,0]*dx, velocity.v


