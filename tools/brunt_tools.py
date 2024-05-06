
from starter2 import *
import dtools.volavg as volavg
reload(ps)

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

def sigmas_2donly(self):
    #this works.  Not normalized, though.
    self.sigma_x2d = np.sqrt(((self.rho2)**2).sum().real)
    self.sigma_k2d = np.sqrt(self.ps2.power[1:].sum().real)
    self.sigma_k2dk= np.sqrt(( self.ps2.kcen*self.ps2.power)[1:].sum())
    self.Rinv = self.sigma_k2dk.real/self.sigma_k2d.real
    self.sigma_Brunt = self.sigma_x2d.real*self.Rinv

def plot_brunt(ftool,outname, fitrange=None, method='full', ax=None):

    if method=='full':
        #ftool.sigmas_full()
        sigmas_full(ftool)
    elif method=='norm':
        sigmas_norm(ftool)
    elif method=='range':
        sigmas_range(ftool, fitrange)
    savefig=False
    if ax is None:
        savefig=True
        fig,ax=plt.subplots(1,1)
    if fitrange is None:
        mask = slice(None)
    else:
        mask = slice(fitrange[0],fitrange[1])
    ax.plot(ftool.ps2.kcen,      ftool.ps2.power,c=[0.5]*3, label='P2d')
    ax.plot(ftool.ps2.kcen[mask],ftool.ps2.power[mask],c='r', label='P2d')
    ax.plot(ftool.ps3.kcen[mask],ftool.ps3.power[mask],c='g', label='P3d')
    ax.plot(ftool.ps2.kcen[mask],ftool.ps2.kcen[mask]*ftool.ps2.power[mask],c='b', label = 'k P2d')
    ax.set(xscale='log',yscale='log')
    ax.legend(loc=1)
    error = 1-ftool.sigma_Brunt/ftool.sigma_x3d
    text_x=0.05
    text_y=0.35
    dy=0.05
    ax.text(text_x, text_y-1*dy,"%0.2e sigma_x3d"%ftool.sigma_x3d, transform=ax.transAxes)
    ax.text(text_x, text_y-2*dy,"%0.2e sigma_x2d"%ftool.sigma_x2d, transform=ax.transAxes)
    ax.text(text_x, text_y-3*dy,"%0.2e R "%(1./ftool.Rinv), transform=ax.transAxes)
    ax.text(text_x, text_y-4*dy,"%0.2e sigma_B "%ftool.sigma_Brunt,  transform=ax.transAxes)
    ax.text(text_x, text_y-5*dy,"%0.2e  error  "%error,  transform=ax.transAxes)
    ax.text(text_x, text_y-6*dy,"%0.2e  ratio  "%(ftool.sigma_Brunt/ftool.sigma_x3d),  transform=ax.transAxes)

    if savefig:
        fig.savefig(outname)


def sigmas_full(self):
    #this works.  Not normalized, though.
    self.sigma_x3d = np.sqrt((self.rho**2).sum().real)
    self.sigma_k3d = np.sqrt((self.ps3.power).sum().real)
    self.sigma_x2d = np.sqrt(((self.rho2)**2).sum().real)
    self.sigma_k2d = np.sqrt(self.ps2.power[1:].sum().real)
    self.sigma_k2dk= np.sqrt(( self.ps2.kcen*self.ps2.power)[1:].sum())
    self.Rinv = self.sigma_k2dk.real/self.sigma_k2d.real
    self.sigma_Brunt = self.sigma_x2d.real*self.Rinv
    self.R1 =self.sigma_Brunt/self.sigma_x3d
    #just to check that everything works right, do it with the actual 3d power spectrum.
    #R2 should be 1
    self.Rinv_actual = self.ps3.power.sum()/self.ps2.power.sum()
    self.sigma_Brunt_actual = self.sigma_x2d*self.Rinv_actual
    self.R2 = self.sigma_Brunt_actual/self.sigma_x3d

def sigmas_range(self, fitrange):
    mask = slice(fitrange[0],fitrange[1])
    Nz = self.rho.size
    N2d = self.rho2.size
    self.mean_rho=(self.rho).sum()/Nz
    self.mean_column= self.rho2.sum()/N2d
    self.sigma_x3d  =np.sqrt(((self.rho-self.mean_rho)**2).sum().real/Nz)
    self.sigma_k3d  =np.sqrt((self.ps3.power[mask]).sum().real/Nz)
    self.sigma_k2dk =np.sqrt(( self.ps2.kcen*self.ps2.power)[mask].sum()/Nz)
    self.sigma_x2d  =np.sqrt(((self.rho2-self.mean_column)**2).sum().real/N2d)
    self.sigma_k2d  =np.sqrt(self.ps2.power[mask].sum().real/N2d)

    #self.Rinv = ()/((self.power_1d2)[1:].sum()/N2d)
    self.Rinv = (self.sigma_k2dk/self.sigma_k2d).real
    #this should give us the right answer.
    #Rinv_actual = (self.power_1d3[1:].sum()/Nz)/(self.power_1d2[1:].sum()/N2d)
    self.Rinv_actual = (self.sigma_k3d/self.sigma_k2d).real
    self.sigma_Brunt = (self.sigma_x2d*self.Rinv).real
    self.sigma_Brunt_actual = (self.sigma_k2d*self.Rinv_actual).real

    R1 =self.sigma_Brunt/self.sigma_x3d
    R2 = self.sigma_Brunt_actual/self.sigma_x3d




#import yt
#def get_cubes(sim,frame,do_rho_4=False):
#    ds = yt.load("/data/cb1/Projects/P49_EE_BB/%s/DD%04d/data%04d"%(sim, frame, frame))
#    print('get cg')
#    cg = ds.covering_grid(0, [0.0]*3, [512]*3)
#
#    dds = cg.dds
#    rho_full = cg["density"].v
#    rho = volavg.volavg(rho_full,rank=3,refine_by=2)
#    output = rho_full, rho
#    if do_rho_4:
#        rho_4 = volavg.volavg(rho,rank=3,refine_by=2)
#        output = rho_full, rho, rho_4
#
#    return output

def fake_powerlaw(N,alphaT,kslice,Amplitude=1,phase=False, rando=False):
    kI = np.fft.fftfreq(N)
    kx,ky,kz=np.meshgrid(kI,kI,kI)
    rrr = np.sqrt(kx**2+ky**2+kz**2)
    Ahat = np.zeros_like(rrr*1j)
    k = np.unique(np.abs(kx))
    if kslice is not None:
        rmin = k[kslice].min()
        rmax = k[kslice].max()
    else:
        rmin=0
        rmax=rrr.max()
    ok = (rrr>rmin)*(rrr<=rmax)
    print(ok.sum()/ok.size)
    ok_r=ok

    #alphaT is total power. 2 makes a flat power
    alphaV=alphaT-2 #alphaV is volume avg power.
    alphaprime=(alphaT-2)/2
    Ahat[ok]=Amplitude*rrr[ok]**alphaprime
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

