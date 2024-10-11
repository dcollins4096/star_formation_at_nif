
from starter2 import *
import tools.brunt_tools as bt
reload(bt)
import tools.fake_powerlaw as fp
reload(fp)

if 1:
    N = 128
    alphaT=-0.5     #total power slope
    alphaV=alphaT-2 #avg power slope
    kmin=2
    kmax=-2
    Q = fp.fake_powerlaw(N,alphaT,kmin,kmax, phase=True, isovec=[1,1,1],rando=False, seed = 8675309)
    sx=slice(None,None);sy=slice(None,None);sz=slice(None,None)
    #sx=slice(N//4,3*N//4);sy=slice(N//4,3*N//4);sz=slice(N//4,3*N//4)
    #sx=slice(N//4,3*N//4);sy=slice(N//4,None);sz=slice(N//4,3*N//4)
    sy=slice(N//4,3*N//4);sx=slice(N//4,None);sz=slice(N//4,3*N//4)

    Qp = Q[sx,sy,sz]
    ftool=bt.fft_tool(Qp)
    ftool.do3()
    ftool.do2(projax=0, apodize=0)
    ftool_full=bt.fft_tool(Q)
    ftool_full.do3()
    ftool_full.do2(projax=0, apodize=0)

    bt.plot_set(ftool,outname='%s/fake_powerlaw'%plot_dir, other=ftool_full)

