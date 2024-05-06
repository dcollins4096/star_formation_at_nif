from starter1 import *
import tools.davetools as davetools
def simple_phase(field1,field2,log=False,ax=None,nbins=64, weights=None,bins=None ):
    ext1=davetools.extents()
    ext2=davetools.extents()
    ext1(field1)
    ext2(field2)
    if bins is not None:
        bins1=bins[0]
        bins2=bins[1]
    else:
        if log:
            bins1=np.geomspace(ext1.minmax[0],ext1.minmax[1],nbins)
            bins2=np.geomspace(ext2.minmax[0],ext2.minmax[1],nbins)
        else:
            bins1=np.linspace(ext1.minmax[0],ext1.minmax[1],nbins)
            bins2=np.linspace(ext2.minmax[0],ext2.minmax[1],nbins)
    hist,xbins,ybins=np.histogram2d( field1,field2, bins=[bins1,bins2],weights=weights)
    output=helper( hist, xbins, ybins,ax=ax)
    return output
def contour(h_in,xbins_in,ybins_in, ax=None,transpose=False, levels=None,**contour_args):
    #takes the output of np.histogram2d (or any other 2d histogram)
    #xbins is 1 larger than h.size[0].

    xbins = xbins_in+0
    ybins = ybins_in+0
    h=h_in+0

    if transpose:
        h = h.transpose()
        temp = xbins.transpose()
        xbins=ybins.transpose()
        ybins=temp


    if xbins.size > h.shape[0]:
        xbins = 0.5*(xbins[1:] + xbins[:-1])
        ybins = 0.5*(ybins[1:] + ybins[:-1])

    nx = len(xbins) ; ny=len(ybins)
    TheX = np.r_[(ny)*[xbins]].transpose()
    TheY = np.r_[(nx)*[ybins]]
    ax.contour(TheX,TheY,h, levels, **contour_args)
def helper(h_in,xbins_in,ybins_in, cmap_name = 'viridis', zlim=None, ax=None,transpose=False, **pcolormesh_args):
    #takes the output of np.histogram2d (or any other 2d histogram)
    #xbins is 1 larger than h.size[0].

    xbins = xbins_in+0
    ybins = ybins_in+0
    h=h_in+0

    if transpose:
        h = h.transpose()
        temp = xbins.transpose()
        xbins=ybins.transpose()
        ybins=temp


    if xbins.size > h.shape[0]:
        xbins = 0.5*(xbins[1:] + xbins[:-1])
        ybins = 0.5*(ybins[1:] + ybins[:-1])

    nx = len(xbins) ; ny=len(ybins)
    TheX = np.r_[(ny)*[xbins]].transpose()
    TheY = np.r_[(nx)*[ybins]]

    if zlim is None:
        zmin = h[h>0].min()
        zmax = h.max()
    else:
        zmin = zlim[0]
        zmax = zlim[1]
    norm = mpl.colors.Normalize( vmin =zmin, vmax=zmax)
    norm = mpl.colors.LogNorm( vmin =zmin, vmax=zmax)
    cmap = copy.copy(mpl.cm.get_cmap(cmap_name))
    cmap.set_under('w')
    if 'shading' not in pcolormesh_args:
        pcolormesh_args['shading'] = 'nearest'
    pcolormesh_args['norm']=norm
    pcolormesh_args['cmap']=cmap
    if ax is not None:

        ploot=ax.pcolormesh( TheX, TheY, h, **pcolormesh_args)

    output = {'TheX':TheX, 'TheY':TheY, 'norm':norm,'cmap':cmap, 'plot':ploot, 'h':h}
    return output




