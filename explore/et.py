
from starter2 import *
import tools.shot as shot
import regions.subregions as subregions
import tools.equal_probability_binner as epb
reload(epb)

Nbins=4
a = np.arange(7)
#a = np.array([0,1,2,3,4,5,6,8])
fig,ax=plt.subplots(1)
pdf, cen, wid=epb.equal_prob(a,Nbins,ax=ax)
fig.savefig('%s/thing'%(plot_dir))

