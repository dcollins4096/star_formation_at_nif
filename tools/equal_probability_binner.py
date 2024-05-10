from dtools.starter1 import *
def equal_prob(arr,Nbins,ax=None,cuml=False):
    v = copy.copy(arr)
    Ntotal = v.size
    Npoints = np.ceil(Ntotal/(Nbins)).astype('int')
    Remaining = Ntotal % Npoints
    v.sort()
    ind = np.arange(0,Ntotal,Npoints)
    hist = ind[1:]-ind[:-1]
    edges = v[list(ind)]
    wid = edges[1:]-edges[:-1]

    #the last bin always needs special treatment.
    #Would love to find a more elegant solution
    if Remaining==0:
        Remaining=Npoints
    hist = np.concatenate([hist,[Remaining]])
    last_bin_width = v[-1]-v[ind[-1]-1]
    wid = np.concatenate([wid,[last_bin_width]])

    pdf = hist/(wid*Ntotal)
    cen = edges+0.5*wid

    #print('Ntotal, Nbins, Npoints, left',v.size, Nbins, Npoints, Remaining)
    #print('R',Remaining)
    #print('ind',ind)
    #print('hist',hist)
    #print('wid',wid)
    #print('pdf',pdf)
    #print('one?',(pdf*wid).sum())
    #error = (pdf*wid).sum()*v.size-v.size
    if ax is not None:
        do_log=False
        if pdf.max()/pdf.min() > 1e3:
            do_log=True
        ax.bar(cen,pdf,width=wid)
        if do_log: ax.set(yscale='log')
        if cuml==True:
            y = np.arange(v.size)/v.size
            ax2 = ax.twinx()
            ax2.plot(v,y,c='k')
            if do_log: ax2.set(yscale='log')
    return pdf, cen, wid
