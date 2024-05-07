#
# shot.shot is one frame.  Contains the raw image and converted density.
# shot.device is two frames, and has tools for the full analysis.
#

from starter2 import *
import regions.subregions as subregions
reload(subregions)
import regions.physical_values as phys
reload(phys)
reload(horz)
import tools.brunt_tools as bt
reload(bt)


class shot():
    def __init__(self,name, lines=[200,400], lightmodel=0, smooth=0):
        print('model, shot',lightmodel)
        self.name=name
        self.image, self.rho, self.I0, self.zero_value= phys.image_to_density(name, lightmodel)
        if smooth>0:
            units=self.rho.units
            self.rho = gaussian_filter(self.rho,smooth)
            self.rho =self.rho*units

        sl = slice(*lines)
        self.rho_cut = self.rho[sl,:]
        self.std =    self.rho_cut.std(axis=0)
        self.rhobar = self.rho_cut.mean(axis=0)
        self.x = phys.get_x(name)
        self.lines=lines

class device():
    def __init__(self,name,lightmodel=0, lines=[200,400],smooth=0):
        self.lightmodel=lightmodel
        self.name = name
        self.name1='%s_t1'%name
        self.name2='%s_t2'%name
        self.shot1 = shot(self.name1,lines,lightmodel,smooth)
        self.shot2 = shot(self.name2,lines,lightmodel,smooth)
        self.lines = lines

    def image_density1(self,fname):
        fig,axes=plt.subplots(3,1, figsize=(4,8))
        ax0=axes[0];ax1=axes[1];ax2=axes[2]
        #ax3=axes[0][1];ax4=axes[1][1];ax5=axes[2][1]

        ax0.plot( self.shot1.rho_cut.transpose(), c=[0.5,0.5,0.5,0.1], linewidth=0.1)
        ax0.plot( self.shot2.rho_cut.transpose(), c=[0.5,0.5,0.5,0.1], linewidth=0.1)
        ax0.plot( self.shot1.rhobar,c='r')
        ax0.plot( self.shot2.rhobar,c='b')
        if 1:#
            cmap = copy.copy(mpl.cm.get_cmap('viridis'))
            cmap.set_under('w')
            mmax = self.shot1.rho.max()
            mmin = max([0,self.shot1.rho.min()])
            norm = mpl.colors.Normalize(vmin=mmin,vmax=mmax)
        ax1.imshow( self.shot1.rho,cmap=cmap,norm=norm)
        ax2.imshow( self.shot2.rho,cmap=cmap,norm=norm)
        ax1.axhline( self.shot1.lines[0],c='r')
        ax1.axhline( self.shot1.lines[1],c='r')
        ax2.axhline( self.shot1.lines[0],c='r')
        ax2.axhline( self.shot1.lines[1],c='r')
        ax0.set(xlabel='x [pixel]',ylabel=r'$\rho(x)$')
        ax1.set(xlabel='x [pixel]',ylabel='y [pixel]')
        ax2.set(xlabel='x [pixel]', ylabel='y [pixel]')
        fig.tight_layout()
        fig.savefig('%s/density_1_%s.pdf'%(plot_dir,self.name))
        plt.close(fig)

    def image_density2(self,fname):
        fig,axes=plt.subplots(3,2, figsize=(12,12))
        ax0=axes[0][0];ax1=axes[1][0];ax2=axes[2][0]
        ax3=axes[0][1];ax4=axes[1][1];ax5=axes[2][1]

        ax0.plot( self.shot1.rho_cut.transpose(), c=[0.5,0.5,0.5,0.1], linewidth=0.1)
        ax0.plot( self.shot2.rho_cut.transpose(), c=[0.5,0.5,0.5,0.1], linewidth=0.1)
        ax0.plot( self.shot1.rhobar,c='r')
        ax0.plot( self.shot2.rhobar,c='b')
        if 1:#
            cmap = copy.copy(mpl.cm.get_cmap('viridis'))
            cmap.set_under('w')
            mmax = self.shot1.rho.max()
            mmin = max([0,self.shot1.rho.min()])
            norm = mpl.colors.Normalize(vmin=mmin,vmax=mmax)
        ax1.imshow( self.shot1.rho,cmap=cmap,norm=norm)
        ax2.imshow( self.shot2.rho,cmap=cmap,norm=norm)
        ax4.imshow( self.shot1.image)#,cmap=cmap,norm=norm)
        ax5.imshow( self.shot2.image)#,cmap=cmap,norm=norm)
        ax1.axhline( self.shot1.lines[0],c='r')
        ax1.axhline( self.shot1.lines[1],c='r')
        ax2.axhline( self.shot1.lines[0],c='r')
        ax2.axhline( self.shot1.lines[1],c='r')
        fig.tight_layout()
        fig.savefig('%s/density_%s.pdf'%(plot_dir,self.name))
        plt.close(fig)

    def bumper(self, rng,fix_shift_x=None, fname=None):
        from scipy.interpolate import CubicSpline
        from scipy.optimize import curve_fit
        sl = slice(*rng)
        #take the units off for the fitter.
        x_units = self.shot1.x[sl].units
        rho_units = self.shot1.rhobar.units
        x_hold = self.shot1.x[sl].v
        y_hold = self.shot1.rhobar[sl].v
        x_move = self.shot2.x[sl].v
        y_move = self.shot2.rhobar[sl].v
        interpolator = CubicSpline(x_move,y_move)
        def test_func(x,dx,dy):
            if fix_shift_x is not None:
                my_dx = fix_shift_x.v
            else:
                my_dx = dx
            return interpolator(x+my_dx) + dy
        self.test_func=test_func
        popt,pcov = curve_fit(test_func,x_hold,y_hold,p0=[0,0])

        if fix_shift_x is not None:
            self.shift_x = fix_shift_x
        else:
            self.shift_x = popt[0]*x_units
        self.shift_rho = popt[1]*rho_units
        self.rhobar_2 = np.interp( self.shot2.x.v+popt[0], self.shot2.x.v, self.shot2.rhobar) + popt[1]
        self.rhobar_2 *= rho_units
        self.x_2 = self.shot2.x #+ self.shift_x*x_units #don't shift twice.
        self.rhobar_1 = self.shot1.rhobar
        self.x_1 = self.shot1.x

        if fname is not None:
            fig,axes=plt.subplots(1,2, figsize=(6,3))
            ax0=axes[1];ax1=axes[0]
            off_blue = [0.5,0.5,1,0.5]
            one_zone = x_hold[1]-x_hold[2]
            print('one zone', one_zone, 'shift',popt[0],'shift in zones', popt[0]/one_zone)
            if 1:
                LW = 0.2
                #ax0.plot(y_hold, c='r')
                #ax0.plot(y_move, c=off_blue)
                #ax0.plot(test_func(x_hold,popt[0],popt[1]), c='b')
                ax0.plot( self.shot1.rhobar[sl], c='r',linewidth=LW)
                ax0.plot( self.shot2.rhobar[sl], c=off_blue,linewidth=LW)
                ax0.plot( self.rhobar_2[sl],c='b',linewidth=LW)
                ax1.plot(self.shot1.rhobar, c='r', linewidth=LW)
                ax1.plot(self.shot2.rhobar, c=off_blue, linewidth=LW)
                ax1.plot(self.rhobar_2, c='b', linewidth=LW)
                ax1.axvline(rng[0],c=[0.5]*4, linewidth=LW)
                ax1.axvline(rng[1],c=[0.5]*4, linewidth=LW)
                ax1.set(xlabel=r'$x\ [pixel]$', ylabel=r'$\rho$')
                ax0.set(xlabel=r'$x\ [pixel]$', ylabel=r'$\rho$')
            if 0:
                ax0.plot(self.shot1.x[sl], self.shot1.rhobar[sl], c='r')
                ax0.plot(self.shot2.x[sl], self.shot2.rhobar[sl], c=off_blue)
                ax0.plot(self.x_2[sl], self.rhobar_2[sl],c='b')
                #ax0.plot(x_hold,y_hold, c='r')
                #ax0.plot(x_move,y_move, c=off_blue)
                #ax0.plot( x_hold, test_func(x_hold,popt[0],popt[1]), c='b')
                ax1.plot(self.shot1.x,self.shot1.rhobar, c='r')
                ax1.plot(self.shot2.x,self.shot2.rhobar, c=off_blue)
                ax1.plot(self.x_2, self.rhobar_2, c='b')
                ax1.axvline(self.shot1.x[rng[0]],c=[0.5]*4)
                ax1.axvline(self.shot1.x[rng[1]],c=[0.5]*4)
                ax1.set(xlabel=r'$x\ [\mu m]$', ylabel=r'$\rho$')
                ax0.set(xlabel=r'$x\ [\mu m]$', ylabel=r'$\rho$')

            fig.tight_layout()
            fig.savefig('%s/%s'%(plot_dir,fname))
            plt.close(fig)

    def get_velocity(self,rng, fname=None, nbins=16):

        ok = slice(rng[0],rng[1])
        xa = self.x_1[ok]
        ya = self.rhobar_1[ok]
        xb = self.x_2[ok]
        yb = self.rhobar_2[ok]
        #horz.try2(ya,yb,method=1,fname='t1')
        #I1 = nar(horz.ho(yb,ya))
        #dx1 = I1[:,1]-I1[:,0]
        #vel1 = phys.pixel_to_velocity(dx1)
        I2 = horz.ho2(ya=yb,yb=ya) #apologies for this looking backwards)
        dx2 = I2[:,1]-I2[:,0]
        self.vel_dist = phys.pixel_to_velocity(dx2)
        hist, cen = epb.equal_prob( self.vel_dist, nbins)
        max_ind= np.argmax(hist)
        self.vel = cen[ max_ind]
        a = yb; b=ya
        if fname is not None:
            #horz.try2(yb,ya,method=2,fname=fname)
            fig,axes=plt.subplots(1,3, figsize=(12,4))
            ax0=axes[1];ax1=axes[2];ax2=axes[0]
            ax2.plot(self.rhobar_1, c='r')
            ax2.plot(self.rhobar_2, c='b')
            ax2.axvline(rng[0],c=[0.5]*4)
            ax2.axvline(rng[1],c=[0.5]*4)
            ax0.plot(range(len(a)), a, color="b")
            ax0.plot(range(len(b)), b, color="r")
            dx=[]
            for ii, x, y in I2:
                i=int(ii)
                xs = [i, x]
                dx.append(x-i)
                ys = [a[i], y]
                ax0.plot(xs, ys, "g--")
                #plt.plot(x, y, "r+")
            epb.equal_prob(nar(self.vel_dist), 16, ax=ax1)
            ax1.scatter( self.vel, hist[max_ind], marker='*',color='orange')
            ax1.text( 0.5,0.75, r'$v = %0.1f km/s$'%self.vel,transform=ax1.transAxes)
            ax2.set(xlabel='x [pixel]', ylabel=r'$\rho(x)$')
            ax0.set(xlabel='x [pixel]', ylabel=r'$\rho(x)$')
            ax1.set(xlabel='v [km/s]', ylabel=r'$P(v)$')
            fig.tight_layout()
            fig.savefig('%s/%s'%(plot_dir,fname))
            plt.close(fig)


    def get_sigma_rho(self,rng=[550,720], fname = None):
        sl = slice(rng[0],rng[1])
        post_shock = self.shot2.rho_cut[:,sl]
        ftool = bt.fft_tool(post_shock)
        ftool.rho2 = post_shock
        ftool.do2()
        bt.sigmas_2donly(ftool)
        self.sigma_B = ftool.sigma_Brunt

        print('fname',fname)

        if fname is not None:
            fig,axes=plt.subplots(1,4, figsize=(13,3))
            #ax0=axes[0][0];ax2=axes[0][1]#;ax2=axes[0][2]
            #ax1=axes[1][0];ax3=axes[1][1]#;ax5=axes[0][2]
            ax0=axes[0];ax1=axes[1];ax2=axes[2];ax3=axes[3]

            ax0.plot( self.rhobar_1,c='r')
            ax0.plot( self.rhobar_2,c='b')
            ax0.axvline(rng[0],c='r')
            ax0.axvline(rng[1],c='r')
            ax1.imshow( self.shot2.rho)

            y1,y2 = self.shot2.lines
            x1,x2 = rng
            ax1.plot( [x1,x2,x2,x1,x1], [y1,y1,y2,y2,y1],c='r')

            ax2.imshow(post_shock)
            ax3.plot(ftool.ps2.kcen,ftool.ps2.power)
            ax3.set(xscale='log',yscale='log')
            ax3.text(0.5,0.75,r'$\sigma_B = %0.2e$'%self.sigma_B, transform=ax3.transAxes)




            fig.tight_layout()
            fig.savefig('%s/density_variance_%s.pdf'%(plot_dir,self.name))
            plt.close(fig)


    def get_csound(self, mean_density, fname=None, gamma=5./3):

        sl = slice(mean_density[0],mean_density[1])
        rhosl = self.rhobar_1[sl]
        index = np.argmin(rhosl)
        mean_rho = self.rhobar_1[sl][index]
        sl2 = slice( mean_density[0]+index, mean_density[2])
        peak_rho = np.max( self.rhobar_1[sl2])
        peak_rho_index = np.argmax( self.rhobar_1[sl2]) + mean_density[0]+index

        R = peak_rho/mean_rho
        print(R)
        self.cs = np.sqrt( gamma*self.vel**2*(1/R)*(1-1/R))
        self.mean_rho = mean_rho
        self.peak_rho = peak_rho
        self.R = R

        if fname is not None:
            if 1:
                fig,ax=plt.subplots(1,1,figsize=(4,4))
                ax0=ax
                ax0.plot(self.rhobar_1,c='k')
                ax0.axvline(mean_density[0], c=[0.5]*4)
                ax0.axvline(mean_density[1], c=[0.5]*4)
                ax0.axvline(mean_density[2], c=[0.5]*4)

                ax0.scatter(index+mean_density[0], mean_rho, c='r',marker='*',s=100)
                ax0.scatter(peak_rho_index, peak_rho, c='r', marker='*',s=100)

                ax0.text(0.1,0.8,  r'$\frac{\rho_{max}}{\rho_{min}}=%0.2e$'%R, transform=ax0.transAxes)
                ax0.text(0.1,0.75,  r'$c_s=%0.2f km/s$'%self.cs, transform=ax0.transAxes)


                fig.tight_layout()
                fig.savefig('%s/%s'%(plot_dir,fname))

            if 0:
                #long plot.
                fig,axes=plt.subplots(1,3)
                ax0=axes[0]; ax1=axes[1]; ax2=axes[2]
                ax0.plot(self.rhobar_1,c='k')
                ax0.axvline(mean_density[0], c=[0.5]*4)
                ax0.axvline(mean_density[1], c=[0.5]*4)
                ax0.axvline(mean_density[2], c=[0.5]*4)
                ax1.plot(  self.rhobar_1[sl])
                ax1.axvline(index, c=[0.5]*4)
                ax1.axhline(mean_rho, c=[0.5]*4)
                ax2.plot( self.rhobar_1[sl2])
                ax2.text(0.1,0.8,  r'$\frac{\rho_{min}}{\rho_{max}}=%0.2e$'%R, transform=ax2.transAxes)
                ax2.text(0.1,0.75,  r'$c_s=%0.2f km/s$'%self.cs, transform=ax2.transAxes)
                ax0.set(xlabel='x [pixel]', ylabel = r'$\rho$')
                ax1.set(xlabel='x [pixel]', ylabel = r'$\rho$')
                ax2.set(xlabel='x [pixel]', ylabel = r'$\rho$')
                fig.tight_layout()
                fig.savefig('%s/%s'%(plot_dir,fname))
    def get_atwood(self):
        self.atwood_number = self.sigma_B/(self.mean_rho+self.sigma_B)
    def get_sigma_v(self):
        if self.name == 'r60' or self.name == 'r120':
            size = {'r60':60, 'r120':120}[self.name]
            self.sigma_v = 0.6*(size/45)**(1./3)*self.atwood_number*self.vel
        else:
            self.sigma_v = 0.6*self.atwood_number*self.vel

