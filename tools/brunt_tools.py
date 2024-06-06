
from starter2 import *
import yt
from dtools.math import volavg
import dtools.davetools as dt
import dtools.math.power_spectrum as ps
reload(ps)


def get_cubes(sim,frame,do_rho_4=False):
    ds = yt.load("/data/cb1/Projects/P49_EE_BB/%s/DD%04d/data%04d"%(sim, frame, frame))
    print('get cg')
    cg = ds.covering_grid(0, [0.0]*3, [512]*3)

    dds = cg.dds
    rho_full = cg["density"].v
    rho = volavg.volavg(rho_full,rank=3,refine_by=2)
    output = rho_full, rho
    if do_rho_4:
        rho_4 = volavg.volavg(rho,rank=3,refine_by=2)
        output = rho_full, rho, rho_4

    return output

def fake_powerlaw(N,alphaT,kmin,kmax,Amplitude=1,phase=False, rando=False):
    kI = np.fft.fftfreq(N)
    kx,ky,kz=np.meshgrid(kI,kI,kI)
    rrr = np.sqrt(kx**2+ky**2+kz**2)
    Ahat = np.zeros_like(rrr*1j)
    k = np.unique(np.abs(kx))

    rmin = k[kmin]
    rmax = k[kmax]

    ok = (rrr>rmin)*(rrr<=rmax)
    ok_r=ok

    #alphaT is total power. 2 makes a flat power
    alphaV=alphaT-2 #alphaV is volume avg power.
    alphaprime=(alphaT-2)/2
    Ahat[ok]=Amplitude*rrr[ok]**alphaprime
    print("Alpha prime",alphaprime)
    #Ahat -= Ahat.min()-4
    if rando:
        rando=np.random.random(Ahat.size)
        rando.shape=Ahat.shape
        Ahat *= rando
    Ahat[0,0,0] =100*Ahat.max()
    #Ahat=np.ones([N,N,N])
    if phase:
        phi = (np.random.random(Ahat.size)-0.5)*np.pi*2
        phi.shape = Ahat.shape
        print('PHI',phi.min(),phi.max())
        Ahatmag = np.abs(Ahat)
        Ahat = Ahatmag*np.cos(phi)+Ahatmag*np.sin(phi)*1j
        Ahat[0,0,0] = np.abs(Ahat[0,0,0])
    H = Ahat.shape[0]//2
    Atotal=Ahat
    Ahat = Ahat[:,:,:H+1]
    Q = np.fft.irfftn(Ahat)
    return Q


def plot_set(ftool,outname):
    fig,ax=plt.subplots(1,2)
    ax0=ax[0];ax1=ax[1]
    ax0.imshow(ftool.rho2,interpolation='nearest',origin='lower')
    plot_brunt(ftool,ax=ax1)
    fig.savefig(outname)
def plot_brunt(ftool,outname=None, fitrange=None, method=None, ax=None):

    sigmas_full(ftool)
    savefig=False
    if ax is None:
        savefig=True
        fig,ax=plt.subplots(1,1)

    if fitrange is None:
        mask = slice(None)
    else:
        mask = slice(fitrange[0],fitrange[1])

    
    M1 = ftool.ps2.power>1e-16
    M2 = M1
    M3 = ftool.ps3.power>1e-16

    ax.plot( ftool.ps2.kcen[M1], ftool.ps2.avgpower[M1],c='m',label='P2d/V2d')
    ax.plot( ftool.ps3.kcen[M1], ftool.ps3.avgpower[M1],c='r',label='P3d/V3d')
    ax.set(xscale='log',yscale='log')
    ax.legend(loc=1)
    error = 1-ftool.ratio_1
    text_x=0.05
    text_y=0.45
    dy=0.07
    ax.text(text_x, text_y-1*dy,"%0.2e sigma_x3d"%ftool.sigma_x3d, transform=ax.transAxes)
    ax.text(text_x, text_y-2*dy,"%0.2e sigma_x2d"%ftool.sigma_x2d, transform=ax.transAxes)
    ax.text(text_x, text_y-3*dy,"%0.2e R "%(1./ftool.Rinv), transform=ax.transAxes)
    ax.text(text_x, text_y-4*dy,"%0.2e sigma_B "%ftool.sigma_Brunt,  transform=ax.transAxes)
    ax.text(text_x, text_y-5*dy,"%0.2e  error  "%error,  transform=ax.transAxes)
    ax.text(text_x, text_y-6*dy,"%0.2e  ratio  "%(ftool.ratio_1),  transform=ax.transAxes)

    if savefig:
        fig.savefig(outname)

def sigmas_2donly(self):
    #this works.  Not normalized, though.
    self.sigma_x2d = np.sqrt(((self.rho2)**2).sum().real)
    self.sigma_k2d = np.sqrt(self.ps2.power[1:].sum().real)
    self.sigma_k2dk= np.sqrt(( self.ps2.kcen*self.ps2.power)[1:].sum())
    self.Rinv = self.sigma_k2dk.real/self.sigma_k2d.real
    self.sigma_Brunt = self.sigma_x2d.real*self.Rinv

def sigmas_full(self):
    #this works.  Probably.
    self.sigma_x3d = np.sqrt((self.rho**2).sum().real)
    self.sigma_k3d = np.sqrt((self.ps3.power).sum().real)
    self.sigma_x2d = np.sqrt(((self.rho2)**2).sum().real)
    self.sigma_k2d = np.sqrt(self.ps2.power.sum().real)
    self.sigma_k2dk= np.sqrt((2* self.ps2.kcen*self.ps2.power).sum())
    self.Rinv = self.sigma_k2dk.real/self.sigma_k2d.real
    self.sigma_Brunt = self.sigma_x2d.real*self.Rinv
    self.ratio_1 =self.sigma_Brunt/self.sigma_x3d
    #just to check that everything works right, do it with the actual 3d power spectrum.
    #R2 should be 1
    self.Rinv_actual = np.sqrt(self.ps3.power.sum()/self.ps2.power.sum())
    self.sigma_Brunt_actual = self.sigma_x2d*self.Rinv_actual
    self.ratio_2 = self.sigma_Brunt_actual/self.sigma_x3d



class fft_tool():
    def __init__(self,rho):
        self.rho=rho
        self.rho2=None
        self.rho2p=None
        self.done2=False
        self.done3=False

    def do3(self):
        self.ps3 = ps.powerspectrum(self.rho)

    def do2(self,projax=0):
        if self.rho2 is None:
            self.rho2=self.rho.sum(axis=projax)
        self.ps2 = ps.powerspectrum(self.rho2)

    def apodize1(self,projax=0):
        self.rho2=self.rho.sum(axis=projax)
        self.rho2/=self.rho2.shape[projax]
        shape=np.array(self.rho2.shape)
        baseshape=(1.*shape).astype('int')
        base = np.zeros(baseshape)
        start = baseshape//2-shape//2

        base[start[0]:(start[0]+shape[0]), start[1]:(start[1]+shape[1])] = self.rho2

        if 0:
            #things that don't quite work
            self.rho2 = base
            from scipy.signal import general_gaussian
            from scipy.signal import convolve2d
            window = np.outer(general_gaussian(baseshape[0],6,3),general_gaussian(baseshape[0],6,3))
            #np.roll(window,baseshape[0]//2,axis=0)
            #np.roll(window,baseshape[0]//2,axis=1)
        if 0:
            #actually work ok
            window = np.zeros(baseshape)
            window[0:3,0:3]=1
        if 1:
            #works pretty well.
            x = np.arange(baseshape[0])
            sigma_conv=2
            g = np.exp(-x**2/(2*sigma_conv**2))**6
            window = np.outer(g,g)

        window/=window.sum()
        self.window=window
        

        #self.rho2 = scipy.convolve(base, window)
        #self.rho2 = convolve2d(base, window)
        #convolve the window function with the base
        a = np.fft.fftn(base)
        b = np.fft.fftn(window)
        c = a*b
        self.rho2 = np.fft.ifftn(c)
        q=np.abs(self.rho2.imag).sum()/self.rho2.imag.size
        if q>1e-13:
            print("!!!!!!!!!!!!!!!!!Imaginary",q)
        self.rho2=self.rho2.real

