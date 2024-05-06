from starter2 import *
import regions.subregions as subregions
reload(subregions)
unit_kappa = unyt.cm**2/unyt.g
unit_density = unyt.g/unyt.cm**3
alpha_B = 0.777*unit_kappa #cm^2/g
alpha_U = 4.697*unit_kappa 
alpha_D = 46.087*unit_kappa 
L_B = 0.06*unyt.cm #cm, 2 layers of 0.03
L_U = 0.14*unyt.cm
L_D = 0.06*unyt.cm
L_F = L_U+L_D
rho_B = 1.845*(unit_density) #g/cm^3
rho_U = 0.090*(unit_density)
rho_D = 0.099*(unit_density)
r = rho_D/rho_U
tau_B = rho_B*L_B*alpha_B
tau_U = rho_U*L_U*alpha_U
tau_D = rho_D*L_D*alpha_D
kappa_bar = 1/(alpha_D*L_D*r+alpha_U*L_U)
t2 = np.exp(alpha_B*L_B*rho_B)

dx_pixel=2.5e-6*unyt.m #micron
t2 = 40e-9*unyt.s #ns
t1 = 38e-9*unyt.s #ns
delta_t = t2-t1


def get_x(shot):
    image = subregions.TNA[shot]
    x = np.arange(image.shape[1])*dx_pixel
    return x
def pixel_to_velocity(dx,dt=(t2-t1),pixel=dx_pixel):
    return (dx*pixel/dt).in_units('km/s')
def compute_I0(I_pixel):
    I_0 = I_pixel * np.exp(tau_B+tau_U+tau_D)
    return I_0
def compute_rho(I_pixel, I_0, zero=0):
    fix1 = (I_pixel-zero)/I_0

    rho = (-np.log(fix1) -tau_B)*kappa_bar
    return rho

def image_to_density(shot, model=0):
    image = subregions.TNA[shot]
    ps = subregions.preshock_region[shot]
    if model == 0:
        I0 = compute_I0(ps).max()
        image_sort = copy.copy(image.flatten())
        image_sort.sort()
        delta = image_sort[1]-image_sort[0] 
        zero_value = image_sort[0]-delta
    elif model == 1:
        I0 = compute_I0(ps).max()
        zero_value = subregions.get_zero(shot)
    elif model == 2:
        I0 = compute_I0(ps).max()
        image_sort = copy.copy(image.flatten())
        image_sort.sort()
        delta = image_sort[1]-image_sort[0] 
        zero_value = image_sort[0]-delta
        #zero_value = 2*subregions.get_zero(shot)
    Q = compute_rho(image, I0, zero=zero_value)
    return image, Q

