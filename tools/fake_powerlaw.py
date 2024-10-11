
from starter2 import *
def fake_powerlaw(N,alphaT,kmin,kmax,Amplitude=1,phase=False, rando=False, isovec=[1,1,1],seed=None):
    kI = np.fft.fftfreq(N)
    kx,ky,kz=np.meshgrid(kI,kI,kI, indexing='ij')
    rrr = np.sqrt(isovec[0]*kx**2+isovec[1]*ky**2+isovec[2]*kz**2)
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
    if seed:
        print('hey!')
        np.random.seed(seed)
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

