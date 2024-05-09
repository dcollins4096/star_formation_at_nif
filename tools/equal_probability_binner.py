from dtools.starter1 import *
def equal_prob(arr,Nbins,ax=None,cuml=False):
    v = copy.copy(arr)
    Ntotal = v.size
    Npoints = np.ceil(Ntotal/(Nbins)).astype('int')
    Remaining = Ntotal % Nbins
    v.sort()
    ind = np.arange(0,Ntotal,Npoints)
    hist = ind[1:]-ind[:-1]
    edges = v[list(ind)]
    wid = edges[1:]-edges[:-1]

    #the last bin needs special handling
    hist = np.concatenate([hist,[Remaining]])
    last_bin_width = v[-1]-v[ind[-1]-1]
    wid = np.concatenate([wid,[last_bin_width]])

    pdf = hist/(wid*Ntotal)
    cen = edges+0.5*wid

    #print('Ntotal, Nbins, Npoints, left',v.size, Nbins, Npoints, Remaining)
    #print('hist',hist)
    #print('wid',wid)
    #print('pdf',pdf)
    #print('one?',(pdf*wid).sum())
    #error = (pdf*wid).sum()*v.size-v.size
    if ax is not None:
        ax.bar(cen,pdf,width=wid)
        if cuml==True:
            y = np.arange(v.size)/v.size
            ax2 = ax.twinx()
            ax2.plot(v,y,c='k')
    return hist, cen
