
from starter2 import *
import tools.shot as shot
import regions.subregions as subregions
import tools.equal_probability_binner as epb
reload(shot)
reload(epb)

name='s120_t2'
sh = shot.shot(name,lightmodel=3,smooth=3) 
fig,ax=plt.subplots(2,2)
ax0=ax[0][0];ax1=ax[0][1]
ax2=ax[1][0];ax3=ax[1][1]

cmap = copy.copy(mpl.colormaps.get_cmap('viridis'))
cmap.set_over('w')

if 1:
    norm = dt.norm_extrema(sh.image.flatten())
    image=sh.image
    image = gaussian_filter(image,10)
    pre=subregions.preshock_region[name]

    norm1 = mpl.colors.Normalize(vmin=0,vmax=400)
    norm1 = dt.norm_extrema(subregions.preshock_region[name].flatten(),frac=0.05)
    #norm1 = mpl.colors.Normalize(vmin=0,vmax=200)
    p=ax0.imshow(image,norm=norm1,cmap=cmap)
    fig.colorbar(p,ax=ax0)

    p1=ax1.imshow(pre,norm=norm1,cmap=cmap)
    pdf,cen,wid=epb.equal_prob(image.flatten(),32,ax=ax2,cuml=-1)
    pdf2,cen2,wid2=epb.equal_prob(pre.flatten(),16,ax=ax3)
    fig.colorbar(p1,ax=ax1)
    fig.savefig('%s/norm_%s'%(plot_dir,name))
    plt.close(fig)
