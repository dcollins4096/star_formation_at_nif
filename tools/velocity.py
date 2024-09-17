from dtools.starter1 import *
import tools.horizontal_distance as horz
reload(horz)


def vel1(rho1, rho2, dx=1, dt=1, ax=None):
    """get the velocity between two densities using lines at constant density"""
    I2 = horz.ho2(ya=rho1, yb=rho2)
    dx2 = I2[:,1]-I2[:,0]
    velocity = dx2*dx/dt
    if ax is not None:
        thex = np.arange(len(rho1))*dx
        ax.plot(thex, rho1)
        ax.plot(thex, rho2)
        for ii, x, y in I2:
            i=int(ii)
            xs = [i*dx, x*dx]
            ys = [rho1[i], y]
            ax.plot(xs, ys, "g--")
        ax2=ax.twinx()
        ax2.plot(I2[:,0]*dx, velocity)

    return velocity.v


