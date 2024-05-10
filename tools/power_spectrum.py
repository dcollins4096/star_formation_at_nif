from starter2 import *
import scipy.stats

class powerspectrum():
    def __init__(self,arr):
        Nhat=np.fft.fftn(arr)
        rhohat = np.abs(Nhat)**2
        rhohat /= rhohat.size
        kx = np.fft.fftfreq(rhohat.shape[0])
        kabs = np.sort(np.unique(np.abs(kx)))
        rank = len(arr.shape)
        if rank == 2:
            kkx,kky=np.meshgrid(kx,kx)
            k = np.sqrt(kkx**2+kky**2)
        elif rank == 3:
            kkx,kky,kkz=np.meshgrid(kx,kx,kx)
            k = np.sqrt(kkx**2+kky**2+kkz**2)
        power, bins, counts =scipy.stats.binned_statistic(k.flatten(), rhohat.flatten(), bins=kabs,statistic='sum')
        bc = 0.5*(bins[1:]+bins[:-1])
        self.Nhat=Nhat  
        self.rho=arr
        self.rhohat=rhohat
        self.k = k
        self.power=power.real
        self.kcen=bc

#import fourier_tools_py3.fourier_filter as Filter
#class powerspectrum_old():
#    def __init__(self,array):
#        self.array=array
#        self.fft = np.fft.fftn( self.array )
#        self.power=self.fft*np.conjugate(self.fft)
#        self.power/=self.power.size
#        ff = Filter.FourierFilter(self.power)
#        self.power_1d = np.array([self.power[ff.get_shell(bin)].sum() for bin in range(ff.nx)])
#        self.Nzones = np.array([ff.get_shell(bin).sum() for bin in range(ff.nx)])
#        self.kcen=ff.get_shell_k()
