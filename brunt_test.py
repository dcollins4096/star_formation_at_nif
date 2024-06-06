
from starter2 import *
import tools.brunt_tools as bt
reload(bt)

if 1:
    N = 128
    alphaT=-1.5     #total power slope
    alphaV=alphaT-2 #avg power slope
    kmin=2
    kmax=-2
    Q = bt.fake_powerlaw(N,alphaT,kmin,kmax, phase=True)
    ftool=bt.fft_tool(Q)
    ftool.do3()
    ftool.do2(projax=0)

    bt.plot_set(ftool,outname='%s/fake_powerlaw'%plot_dir)

