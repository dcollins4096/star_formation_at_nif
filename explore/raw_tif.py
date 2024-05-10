
from starter2 import *

import regions.tiff_poker as TP
import regions.regions_coarse as regions_coarse
import regions.regions_align as regions_align
reload(regions_coarse)

#regions_coarse.raw_image(3,'shot_day_2_0')
#regions_coarse.raw_image(4,'shot_day_2_1')
#regions_coarse.raw_image(5,'shot_day_2_2')

V = TP.viewer(which=4)
#V.guess_scale()
AD=V.all_data
fig,ax=plt.subplots(1,1)
norm = mpl.colors.LogNorm(vmin=AD.min(),vmax=AD.max())
p=ax.imshow(AD,norm=norm)
fig.colorbar(p)
fig.savefig('%s/%s'%(plot_dir,'tmp2'))



