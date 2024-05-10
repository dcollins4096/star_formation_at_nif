
from starter2 import *

def image_pre_align(A1,A2,fname):
    fig,ax=plt.subplots(2,1)

    ext = dt.extents()
    ext(A1)
    ext(A2)
    #norm = mpl.colors.Normalize(vmin=ext.minmax[0],vmax=ext.minmax[1])
    norm = dt.norm_extrema(np.concatenate([A1.flatten(),A2.flatten()]))
    ax[0].imshow(A1,norm=norm)
    ax[1].imshow(A2,norm=norm)
    #print(A1.max())
    #print(A2.max())
    #ax[2].hist(A1.flatten(),histtype='step')
    #ax[2].hist(A2.flatten(),histtype='step')
    fig.savefig('%s/%s'%(plot_dir,fname))

def align(A1,A2,name_base='shot', xrange=[290,350], thresh=250, x_override=None, fid_line=600, preserve_x=False):
    #Find the fiducial.
    #Look for the point of the fiducial.
    #Trim based on the tip of the fiducial.
    fig,axes=plt.subplots(4,2, figsize=(12,12))
    ax0=axes[0][0];ax1=axes[0][1]
    ax2=axes[1][0];ax3=axes[1][1]
    ax4=axes[2][0];ax5=axes[2][1]
    ax6=axes[3][0];ax7=axes[3][1]

    Nx1,Ny1=A1.shape
    Nx2,Ny2=A2.shape


    #lines in Y around the fiducial.
    #look for the highest minimum
    fid_edges_1=[]
    fid_x=np.arange(xrange[0],xrange[1])
    plot_locators=True
    for vfid in fid_x:
        line=A1[:,vfid]
        fid_edge = np.where(line>thresh)[0].max()
        #fid_edge = np.where(line>thresh)[0].max()
        fid_edges_1.append(fid_edge)
        if plot_locators:
            ax4.plot(line,linewidth=0.1)
            ax0.axvline(vfid,linewidth=0.1,c='r')
    fid_x_1 = fid_x[ np.argmin(fid_edges_1)]
    fid_y_1 = min(fid_edges_1)

    #repeat for the 2nd shot
    fid_edges_2=[]
    fid_x=np.arange(xrange[0],xrange[1])
    for vfid in fid_x:
        line=A2[:,vfid]
        fid_edge = np.where(line>thresh)[0].max()
        fid_edges_2.append(fid_edge)
        if plot_locators:
            ax5.plot(line,linewidth=0.1)
            ax1.axvline(vfid,linewidth=0.1,c='r')
    fid_x_2 = fid_x[ np.argmin(fid_edges_2)]
    fid_y_2 = min(fid_edges_2)

    ax4.axhline(thresh)
    ax5.axhline(thresh)

    #find the shift
    if x_override is not None:
        dx = x_override
    else:
        dx = (fid_x_2-fid_x_1)
    dy = (fid_y_2-fid_y_1)

    #extents.  We intersect the two images.
    dy1 =  np.abs(min([dx,0]))
    dy2 =  np.abs(min([-dx,0]))
    dx1 =  np.abs(min([dy,0]))
    dx2 =  np.abs(min([-dy,0]))


    LLLx = min([ Nx1-dx1, Nx2-dx2])
    LLLy = min([ Ny1-dy1, Ny2-dy2])
    if preserve_x:
        T1 = A1[dx1:dx1+LLLx, dy1:]
        T2 = A2[dx2:dx2+LLLx, dy2:]
    else:
        T1 = A1[dx1:dx1+LLLx, dy1:dy1+LLLy]
        T2 = A2[dx2:dx2+LLLx, dy2:dy2+LLLy]

    norm = dt.norm_extrema(np.concatenate([A1.flatten(),A2.flatten()]))
    ax0.imshow(A1, origin='lower', interpolation='nearest',norm=norm)
    ax1.imshow(A2, origin='lower', interpolation='nearest',norm=norm)
    hfid=fid_line
    ax0.axhline(hfid)
    ax1.axhline(hfid)

    ax2.plot( A1[hfid,:])
    ax2.plot( A2[hfid,:])
    ax2.set(title='x at %d'%hfid)
    ax3.set(title='blank')

    ax6.plot(fid_edges_1)
    ax7.plot(fid_edges_2)
    ax6.set(title='upper fiducial profile')
    ax7.set(title='upper fiducial profile')
    ax4.set(title='y lines around fiducial')
    ax5.set(title='y lines around fiducial')

    ax0.axvline(fid_x_1, linewidth=0.2,c='k')
    ax0.axhline(fid_y_1, linewidth=0.2,c='k')
    ax1.axvline(fid_x_2, linewidth=0.2,c='r')
    ax1.axhline(fid_y_2, linewidth=0.2,c='r')



    output='%s/%s_align1'%(plot_dir,name_base)
    fig.tight_layout()
    fig.savefig(output)
    print(output)
    plt.close('fig')

    if 0:
        def scaler(img):
            out = img - img.min()
            out /= out.max()
            out = out.transpose()
            out *= 255
            out = out.astype('int')
            return out
        #color channel game.
        Zero = np.zeros_like(T1)
        I1 = scaler(np.stack([T1,0.5*T1,Zero]))
        I2 = scaler(np.stack([Zero,0.5*T2,T2]))
        #ax5.imshow(I1, origin='lower',interpolation='nearest')
        #ax6.imshow(I2, origin='lower',interpolation='nearest')
        fig, ax=plt.subplots(1,1, figsize=(12,12))
        tots = I1+I2
        tots = (tots/tots.max() * 255    ).astype('int')
        ax.imshow(tots,interpolation='nearest')
        fig.savefig('%s/%s_aligned_colortest'%(plot_dir,name_base))

    norm = dt.norm_extrema(np.concatenate([T1.flatten(),T2.flatten()]))
    name1 = '%s/%s_shot1_trim_align.png'%(plot_dir,name_base)
    name2 = '%s/%s_shot2_trim_align.png'%(plot_dir,name_base)
    mpl.image.imsave(name1,T1,vmin=norm.vmin, vmax=norm.vmax)
    mpl.image.imsave(name2,T2,vmin=norm.vmin, vmax=norm.vmax)

    return T1, T2





