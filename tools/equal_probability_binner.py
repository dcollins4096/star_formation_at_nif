from dtools.starter1 import *
def equal_prob(arr,Nbins,ax=None):
    v = copy.copy(arr)
    Npoints = v.size//Nbins
    v.sort()
    ind = np.arange(0,v.size,Npoints)
    #ind = np.concatenate([ind,v.size-1])
    edges = v[list(ind)]
    cen = 0.5*(edges[1:]+edges[:-1])
    wid = edges[1:]-edges[:-1]
    hist = 1/wid
    if ax is not None:
        ax.bar(cen,hist,width=wid)
    return hist, cen
