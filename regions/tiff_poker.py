
from starter2 import *

#base_dir="/Users/dcollins/Dropbox/RESEARCH5/Paper68/Data_analysis/Raw_radiographs/play"
#base_dir="/Users/davidcollins/Dropbox/RESEARCH5/Paper68/Data_analysis/Raw_radiographs/play"
base_dir = 'data'
i1="TD_TC090-124_HGXD_IMAGE_N220712-002-999_DROOP_CORR_422421128478532_20220907115837893.tif"
i2="TD_TC090-124_HGXD_IMAGE_N220713-001-999_DROOP_CORR_745405306609760_20220907115905225.tif"
i3="TD_TC090-124_HGXD_IMAGE_N220714-001-999_DROOP_CORR_766909572834217_20220907115956892.tif"
fnames=[i1,i2,i3]

for fn in fnames:
    fullname='%s/%s'%(base_dir,fn)
    if not os.path.exists(fullname):
        print("Error: file not found: ",fullname)

def trimmer(arr,sigma_n=0,fname='imag', vmin=None,vmax=None):
    y = arr.mean(axis=1)
    ok = np.where(y>sigma_n)
    ymin = ok[0].min()
    ymax=ok[0].max()
    x = arr.mean(axis=0)
    ok = np.where(x>sigma_n)
    xmin = ok[0].min()
    xmax=ok[0].max()
    trim = arr[ymin:ymax,xmin:xmax]
    fig,axes=plt.subplots(2,2,figsize=(12,12))
    #ax0=axes[0];ax1=axes[1]
    ax0=axes[0][0];ax1=axes[0][1]
    ax2=axes[1][0];ax3=axes[1][1]
    if vmin is not None:
        vmax=arr.max()
        vmin=0
    norm=mpl.colors.Normalize(vmin=vmin,vmax=vmax)
    ax0.imshow(arr,norm=norm,origin='lower')
    ax1.imshow(trim,norm=norm,origin='lower')
    ax2.plot(y)
    ax2.axhline(sigma_n)
    ax2.axvline(ymin)
    ax2.axvline(ymax)
    ax3.axvline(xmin)
    ax3.axvline(xmax)
    ax3.plot(x)
    ax3.axhline(sigma_n)

    ax0.axhline(ymin)
    ax0.axhline(ymax)
    ax0.axvline(xmin)
    ax0.axvline(xmax)
    fig.savefig(fname)
    plt.close(fig)
    return trim


class viewer():
    def __init__(self,which=None,arr=None):
        if type(which) == int:
            self.fname = "%s/%s"%(base_dir,fnames[which])
            self.all_data = imread(self.fname)
        elif arr is not None:
            self.all_data=arr
        else:
            pdb.set_trace()
            print("Error, not an acceptable choice")
        self.XX = np.arange(self.all_data.shape[1])
        self.YY = np.arange(self.all_data.shape[0])
        self.X1, self.Y1 = np.meshgrid(self.XX,self.YY)

    def image1(self,fname='image0'):

        self.guess_scale()
        vmin,vmax=self.minmax
        norm=mpl.colors.Normalize(vmin=vmin,vmax=vmax)
        fig,axes=plt.subplots(1,1)
        axes.pcolormesh(self.X1,self.Y1,self.all_data,norm=norm)
        axes.set_aspect('equal')
        fig.savefig(fname)
        plt.close(fig)




    def guess_scale(self):
        ad=self.all_data
        test = ad[ad.shape[0]//4:int(3/4*ad.shape[0]),ad.shape[1]//4:int(3/4*ad.shape[1])]
        self.minmax=np.array([test[test>0].min(),test.max()])



    def xtract(self,a=None,b=None,c=None,d=None):
        S1 = slice(c,d)
        S2 = slice(a,b)
        TheX,TheY,TheZ=self.X1[S1,S2],self.Y1[S1,S2],self.all_data[S1,S2]
        return TheX, TheY, TheZ

    def xtract_and_image(self,a=None,b=None,c=None,d=None,vmin=None,vmax=None, fname='image.png',zero=False):
        #ax0=axes[0][0];ax1=axes[0][1]
        #ax2=axes[1][0];ax3=axes[1][1]
        S1 = slice(c,d)
        S2 = slice(a,b)
        TheX,TheY,TheZ=self.X1[S1,S2],self.Y1[S1,S2],self.all_data[S1,S2]
        if vmin is None and vmax is None:
            vmin = TheZ.min()
            vmax = TheZ.max()
        n_plots = 2
        if zero: n_plots=3
        fig,axes=plt.subplots(1,n_plots,figsize=(12,12))
        #for ax in axes.flatten():
        #    ax.set_aspect('equal')
        #ax0=axes
        ax0=axes[0];ax1=axes[1]
        ax0.plot([a,b,b,a,a],[c,c,d,d,c],c='r')
        cmap='viridis'
        if zero:
            maxmax=max([np.abs(vmin),np.abs(vmax)])
            vmax = maxmax
            vmin = -maxmax
            cmap = 'seismic'

        #norm=mpl.colors.Normalize(vmin=vmin,vmax=vmax)
        norm=mpl.colors.Normalize(vmin=vmin,vmax=0)
        ax0.pcolormesh(self.X1,self.Y1,self.all_data,norm=norm)
        p=ax1.pcolormesh(TheX,TheY,TheZ,norm=norm, cmap=cmap)
        fig.colorbar(p,ax=ax1)
        if zero and 1:
            ax2=axes[2]
            if 0:
                bins = np.arange(-50,150,10)
                ax2.hist(TheZ.flatten(), histtype='step', density=False,bins=bins)
                ax2.set(yscale='log')
            if 1:
                #hist,cen=epb.equal_prob(TheZ.flatten(), 16,ax=ax2)
                hist,cen=epb.equal_prob(self.all_data.flatten(), 16,ax=ax2)
                zero = cen[np.argmax(hist)]
                ax2.text(0.5,0.75,zero, transform=ax2.transAxes)

            #dt.phist(TheZ.flatten())
        if zero and 1:
            ax2=axes[2].twinx()
            the_x = TheZ.flatten()+0
            the_x.sort()
            the_y = np.arange(the_x.size)/the_x.size
            ok = np.argmin(np.abs(the_x))
            ax2.axhline(the_y[ok])
            ax2.axvline(0)
            ax2.plot(the_x,the_y)


        fig.tight_layout()
        fig.savefig(fname)
        plt.close(fig)
        return TheX,TheY,TheZ






