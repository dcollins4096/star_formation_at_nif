import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
import types
import glob
import os.path
import tarfile
import h5py
class extents():
    def __init__(self, array=None):
        self.minmax=[]
        self.errors=[]
        if array is not None:
            self(array)

    def __call__(self,array):
        if hasattr(array,'min'):
            the_min=array.min()
            the_max=array.max()
        else:
            the_min=min(array)
            the_max=max(array)
        if len(self.minmax):
            self.minmax[0] = min([the_min,self.minmax[0]])
            self.minmax[1] = max([the_max,self.minmax[1]])
        else:
            self.minmax=[the_min,the_max]
    def __getitem__(self,index):
        return self.minmax[index]
    def __str__(self):
        if len(self.minmax) == 2:
            out= "[%0.2e, %0.2e]"%tuple(self.minmax)
        else:
            out = "undef"
        return out
    def __repr__(self):
        if len(self.minmax) == 2:
            out= "[%0.2e, %0.2e]"%tuple(self.minmax)
        else:
            out = "undef"
        return out
    def check(self,array):
        if array.min() < self.minmax[0]:
            self.errors.append(array.min())
            print("Extent error: min")
        if array.max() > self.minmax[1]:
            self.errors.append(array.max())
            print("Extent error: max")

def read_fft(fname,setname):
    h5ptr = h5py.File(fname,'r')
    output = h5ptr[setname][:]
    return output

def axbonk(ax,xscale='linear',yscale='linear',xlabel='X',ylabel='Y',xlim=None,ylim=None,
          linthreshx=0.1,linthreshy=0.1):
    ax.set_xscale(xscale)
    ax.set_yscale(yscale)
    if xscale == 'symlog':
        ax.set_xscale(xscale,linthreshx=linthreshx)
    if yscale == 'symlog':
        ax.set_yscale(yscale,linthreshy=linthreshy)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
def lim_down(value):
    return 10**(np.floor(np.log10(value)))
def lim_up(value):
    return 10**(np.ceil(np.log10(value)))
def dumb_plt(plot,X,Y,xlabel,ylabel,outname,scale=('linear','linear'), clobber=False, scatter=False,**kwargs):
    if scatter:
        verb = plot.scatter
    else:
        verb= plot.plot
    if clobber:
        plot.clf()

    if X is None:
        X = np.arange(Y.size)
    print("WTF")
    if hasattr(plot,'savefig'):
        output=verb(X,Y,**kwargs)
        plot.xscale(scale[0])
        plot.yscale(scale[1])
        plot.xlabel(xlabel)
        plot.ylabel(ylabel)
        plot.savefig(outname)
    else:
        output=verb(X,Y,**kwargs)
        plot.set_xscale(scale[0])
        plot.set_yscale(scale[1])
        plot.set_xlabel(xlabel)
        plot.set_ylabel(ylabel)
        plot.figure.savefig(outname)
    print(scale)
    plot.yscale('log')
    print(outname)
    return output
    #return line
def collect_extrema(a,b=None):
    """b collects the extrema of a"""
    if b is None:
        b=np.array([min(a),max(a)])
    b[0] = min([min(a),min(b)])
    b[1] = max([max(a),max(b)])
    return b


def ensure_list(obj):
    """
    This function ensures that *obj* is a list.  Typically used to convert a
    string to a list, for instance ensuring the *fields* as an argument is a
    list.
    """
    if obj is None:
        return [obj]
    if type(obj) is not list:
        return [obj]
    return obj

def quarts(array, lower=0.25, upper=0.75, take_log=True):
    print("quarts doesn't work")
    ok = array >0
    pdb.set_trace()
    if take_log:
        hist_for_cbar = np.histogram( np.log10( array[ok].flatten() ),bins=100, normed=True )
    else:
        hist_for_cbar = np.histogram( array[ok].flatten() ,bins=100, normed=True )

    cum = np.cumsum(hist_for_cbar[0])/hist_for_cbar[0].sum()
    first_quart = cum-lower
    if (first_quart == 0).any():
        first_ind = np.where( first_quart == 0 ) [0]
    else:
        first_ind = np.where( first_quart[:-1]*first_quart[1:] <= 0 ) [0]
    last_quart = cum-upper
    if (last_quart == 0).any():
        last_ind = np.where(last_ind == 0 )[0]
    else:
        last_ind = np.where( last_quart[:-1]*last_quart[1:] < 0 ) [0]
    if take_log:
        low_value = 10**np.floor( hist_for_cbar[1][first_ind] )
        high_value = 10**np.ceil(  hist_for_cbar[1][last_ind] )
    else:
        low_value  = np.floor( hist_for_cbar[1][first_ind] )
        high_value = np.ceil(  hist_for_cbar[1][last_ind] )
    return low_value, high_value


def tff(g_code=None,pf=None):
    if pf is not None:
        g_code = pf['GravitationalConstant']

    return np.sqrt(3.*np.pi*4*np.pi/(32*g_code))

def sci(number):
    return "%0.2e"%(number)

def to_tar_gz(source_dir, destination):
    """
    param source_dir: Source directory name.
    param destination: Destination filename.
    (TAR-GZ-Archive *.tar.gz)
    Stolen from http://www.velocityreviews.com/forums/t370204-create-tarfile-using-python.html
    """

    t = tarfile.open(name = destination, mode = 'w:gz')
    t.add(source_dir, os.path.basename(source_dir))
    t.close()

    return True

def no_whites(something):
    """Removes any empty charachters"""
    out = []
    for s in something:
        if s != '':
            out.append(s)
    return out

import shutil
def no_trailing_comments(file):
    nbackups = len(glob.glob('%s.backup*'%file))
    backup_name = '%s.backup%d'%(file,nbackups)
    shutil.move(file, backup_name)
    iptr = open(backup_name,'r')
    optr = open(file,'w')
    for line in iptr:
        if "//" in line:
            optr.write(line[:line.index("//")]+"\n")
        else:
            optr.write(line)
    iptr.close()
    optr.close()

def ImRoot():
    """
    This function accepts a *func*, a set of *args* and *kwargs* and then only
    on the root processor calls the function.  All other processors get "None"
    handed back.
    """
    from yt.config import ytcfg
    if not ytcfg.getboolean("yt","__parallel"):
        return True
    try:
        if ytcfg.getint("yt", "__parallel_rank") > 0: return False
    except:
        if ytcfg.getint("yt","__topcomm_parallel_rank") > 0: return False
    return True

def rainbow_01():
    norm = mpl.colors.Normalize()
    norm.autoscale([0.0,1.0])
    cmap = mpl.cm.jet
    color_map = mpl.cm.ScalarMappable(norm=norm,cmap=cmap)
    return  color_map.to_rgba

class rainbow_map():
    def __init__(self,n):
        norm = mpl.colors.Normalize()
        norm.autoscale(np.arange(n))
        #cmap = mpl.cm.jet
        self.color_map = mpl.cm.ScalarMappable(norm=norm,cmap='jet')
    def __call__(self,val,n_fields=0):
        this_value = self.color_map.to_rgba(val)
        if n_fields > 0:
            this_value = [this_value]*n_fields
        return this_value

    #return  color_map.to_rgba

try:
    algae = mpl.cm.get_cmap("algae")
except:
    pass
def algae_map(n):
    norm = mpl.colors.Normalize()
    norm.autoscale(np.arange(n))
    #cmap = mpl.cm.__dict__['algae']
    color_map = mpl.cm.ScalarMappable(norm=norm,cmap='algae')
    return  color_map.to_rgba

def trans(array,dim):
    """returns the elements of array that aren't dim
    array must be a numpy array.
    trans([a,b,c],1) -> [a,c]"""
    return array[filter(lambda x: x != dim,range(len(array)) ) ]

def relerr(a,b):
    """returns (b-a)/<a,b>"""
    return (b-a)/(0.5*(a+b))

def relerr0(a,b):
    """returns (b-a)/<a,b>"""
    return (b-a)/a

def psave(ax,name):
    print( name)
    ax.savefig(name)

def dsave(ax,name,field_name=None,pf_list=None,tar=True,script_name=None,plt_format='dave'):
    if plt_format != 'dave':
        if plt_format not in name:
            name += "."+plt_format
        psave(ax,name)
        return

    basename = name.split(".")
    #pdb.set_trace()
    if len(basename) > 1:
        format = basename[-1]
        basename = basename[:-1]
        basename = ("%s."*len(basename)%tuple(basename))[:-1] #the last [:-1] is for the extra ."
    else:
        basename = basename[0]
        format = 'whatever'

    if format not in ['png','pdf','eps','whatever','dave']:
        print("Unclear format string:",format)
    directory = "figs/"+basename
    if field_name is not None:
        directory += "_%s"%field_name

    try:
        os.mkdir("figs")
    except:
        pass
    try:
        os.mkdir(directory)
    except:
        pass
    for form in ['png','eps','pdf']:
        if hasattr(ax,'savefig'):
            ax.savefig("%s/%s.%s"%(directory,basename,form))
        if hasattr(ax,'save'):
            ax.save("%s/%s"%(directory,basename), format=form)

    d_html(directory,basename,field_name,pf_list,script_name)
    print("open %s/%s.png"%(directory,basename))
    print("get %s.tar"%(directory))
    to_tar_gz(directory,directory+".tar")


def d_html(directory,basename,field_name,pf_list=None,script_name=None):
    """write meta html."""
    field_to_add = ""
    if field_name is not None:
        field_to_add = "_%s"%field_name
    pname = basename+field_to_add+".png"
    bname = basename
    hname = basename+field_to_add+".html"
    dname = directory
    full_hname = "%s/%s"%(dname,hname)
    print(full_hname)

    hptr = open(full_hname,"w")
    hptr.write("<html><br>\n")
    hptr.write('<!-- In case you want a link:<br>\n')
    hptr.write('<a href = "%s"><br>\n<img src="%s" width=200></a><br>\n--><br>'%(pname,pname))
    hptr.write('<center><img src=%s></center><br>\n'%pname)
    machine = "unknown"
    if 'machine' in os.environ:
        machine = os.environ['machine']
    hptr.write("Plotted on: %s<br>\n"%machine)
    if script_name is not None:
        hptr.write('Script Name:%s<br>\n'%script_name)
    if pf_list is not None and len(pf_list) > 0:
        hptr.write("PF list:<br>\n")
        for pf in pf_list:
            hptr.write("%s<br>"%pf)
    else:
        hptr.write("unknown PF<br>")
    hptr.close()



def dpy(filename,fields):
    """Open hdf5 file *filename* and return *field*, then close the file.
    Collapses 3 function calls into 1."""
    if glob.glob(filename) == []:
        print("No such file", filename)
        return [None]*len(fields)
    fptr = h5py.File(filename,'r')
    output = []
    try:
        if hasattr(fptr,'listnames'):
            field_list = fptr.listnames()
        elif hasattr(fptr,'keys'):
            field_list = fptr.keys()

        for field in ensure_list(fields):
            if field not in field_list:
                print("No field", field, "in file",filename)
                return None
            else:
                if len(fptr[field].shape) > 0:
                    output.append(fptr[field][:])
                else:
                    output.append(fptr[field].value)
    except:
        pdb.set_trace()
        raise
    finally:
        fptr.close()
    return output

def dpy_save(filename,object,fields):
    """Write all *object*[*fields*] to hdf5 file *filename*"""
    fptr = h5py.File(filename,'w')
    for field in fields:
        set = object[field]
        fptr.create_dataset(field,set.shape,data=set)
    fptr.close()

def expform(float_in, format = "%0.1e"):
    """exponent format: 1.5e+03 goes to $1.5 \times 10^{3}$"""
    str1 = format%float_in
    tex_string = "$"
    try:
        exp_pos = str1.index("e")
        exponent_shift = 1
        if str1[exp_pos+1] == "+":
            exponent_shift+1
        mantissa = str1[0:exp_pos]
        exponent = str1[exp_pos+exponent_shift:]
        if float(mantissa) != 1:
            tex_string += mantissa + "\\times "
        tex_string += "10^{" 
        tex_string += str(int(exponent))
        tex_string += "}$"
    except:
        tex_string += str1
        tex_string += "$"
    return '%s'%tex_string

def grep(lookfor,obj):
    if isinstance(obj, types.ListType):
        my_list = obj
    else: my_list = dir(obj)
    for i in my_list:
        if lookfor.upper() in i.upper(): print(i)

def stat(array,strin='', format='%0.16e'):
    template = '['+format+','+format+'] %s %s'
    print(template%(array.min(),array.max(),array.shape,strin))

def nonzerostat(array,strin=''):

    print('[%0.4e,%0.4e] %s %s'%(array[array>0].min(),array[array>0].max(),array.shape,strin))

def morestat(array,strin='',log=False):
    if log:
        mean = 10**(meanRMS(np.log10(array)))
    else:
        mean = meanRMS(array)
    print('[%0.4e,%0.4e \pm %0.4e, %0.4e] %s %s'%(array.min(),mean[0],mean[1],array.max(),array.shape,strin))

def psave(ax,name):
    print(name)
    ax.savefig(name)

def plave(array,filename,axis=None,ax=None,colorbar=True, zlim=None, label="Value",scale='linear',ticks_off=False,
         linthresh=1.):
    """plot and save.  Takes *array*, saves an image and saves it in *filename*
    *zlim* = [a,b] scales the colorbar form a to b."""
    save_fig = False
    if ax is None:
        fig, ax = plt.subplots(1,1)
        save_fig=True

    if len(array.shape) == 2:
        array_2d = array
    else:
        if axis == None:
            print("Must provide an axis: 3d array")
            return
        #array_2d = np.flipud(np.sum( array , axis=axis).transpose())
        array_2d = np.sum( array , axis=axis)

    norm = None
    cmap = mpl.cm.gray
    if zlim is None:
        zlim = [array_2d.min(),array_2d.max()]
    if scale == 'linear':
        norm=mpl.colors.Normalize(vmin=zlim[0], vmax=zlim[1])
    elif scale == 'log':
        norm=mpl.colors.LogNorm(vmin=zlim[0], vmax=zlim[1])
    elif scale == 'symlog':
        norm=mpl.colors.SymLogNorm(linthresh,vmin=zlim[0], vmax=zlim[1])
    plot = ax.imshow(array_2d, interpolation='nearest', origin='lower', norm=norm,cmap=cmap)
    if not ticks_off:
        ax.set_yticks(range(0,array_2d.shape[0],5))
        ax.set_xticks(range(0,array_2d.shape[1],5))
    else:
        ax.set_xticks([])
        ax.set_yticks([])
    print(filename)
    norm = None
    if colorbar:
        colorbar = plt.colorbar(plot, norm=norm)
        colorbar.set_label(label)
        if zlim is not None:
            colorbar.set_clim(zlim[0], zlim[1])


    if save_fig:
        fig.savefig(filename)
        plt.close(fig)
        print(filename)

def wp(axes,number,index_only=False):
    """Which Plot.  Takes 2d list *axes* and returns the row-major *number* element.
    *index_only* returns (i,j)"""
    nx = len(axes[0])
    i = int(number/nx)
    j = number%nx
    if index_only:
        output = i,j
    else:
        output = axes[i][j]
        #output.text(1e-2,1e-1,'bn %d (%d,%d)'%(number,i,j))
    return output

def meanRMS(F):
    M = F.sum()/(F.size)
    rms = np.sqrt(( ( F - M )**2 ).sum()/F.size)
    return nar([M,rms])




def morestat(array,strin='',log=False):
    if log:
        mean = 10**(meanRMS(np.log10(array)))
    else:
        mean = meanRMS(array)
    print('[%0.4e,%0.4e \pm %0.4e, %0.4e] %s %s'%(array.min(),mean[0],mean[1],array.max(),array.shape,strin))

def powerline(this_plt,x1,x2,y1,power,log=True,**kwargs):
    """Plot a powerlaw on the current plot in matplot lib instance *plt*.
    Starts at *x1*, *y1* and runs to *x2*, whatever the power law dictates.
    *log* determines log or linear plot."""
    if False:
        if True:
            x1 = 10**x1
            x2 = 10**x2
        x = [x1,x2]
        b = x1**(-power)*y1
        yf = x2**power*b
        y = [y1, yf]
        if True:
            x1 = np.log10(x1)
            x2 = np.log10(x2)
    if log:
        x = [x1,x2]
        yf = x2**power*x1**(-power)*y1
        y = [y1, yf]
    else:
        x = [x1,x2]
        yf = power*(x2-x1) + y1
        y = [y1,yf]
    this_plt.plot(x,y,**kwargs)

def read_csv(filename):
    file=open(filename,'r')
    lines=file.readlines()
    file.close()
    obs = [ L.split(',') for L in lines]
    values = {}
    for i,n in enumerate(obs[0]):
        values[n] = [ ooo[i] for ooo in obs[1:] ]
        try:
            values[n] = map(float,values[n])
        except:
            pass
        values[n]=nar(values[n])
    return values

def getdata(file_list):
    out = []
    for f in file_list:
        masses = read_csv(f)['Msun']
        out += map(float,masses)                                        
    return nar(out)

def phist(array,width = 8, format = 'f', plot=plt,**kwargs):
    """Accepts an array to  histogram, and prints it in a way that doesn't suck.
    *width* is the width of the output field.
    *format* is the bin format ('f'loat or 'e'xponential)"""
    hist = plot.hist(array,**kwargs)
    val = hist[0]
    bin = hist[1]
    #These are meta-format strings.
    vbase = "%s%dd "%("%",width)
    bbase = "%s%d.2%s "%("%",width,format)
    halfformat = "%s%d.2%s "%("%",width/2,format)
    vformat = ""
    bformat = halfformat
    for n in range(len(val)):
        vformat += vbase
        bformat += bbase
    print(vformat%tuple(val))
    print(bformat%tuple(bin))
    return hist

class ParameterException(Exception):
    def __init__(self, parameter,pf):
        self.value = 'Parameter file %s has no parameter %s'%(str(pf),parameter)
    def __str__(self):
        return repr(self.value)

def tabler(head,rows):
    """Make a latex table."""
    setup = '\\begin{table}[h]' + '\n' + r'\begin{center}' +'\n' +'\caption{}'
    setup += '\\begin{tabular}{ %s }'%('c '*len(head) )
    setup += '\\label{}'+'\n'
    row_format = '%10s'+'& %10s'*(len(head)-1) 
    head = row_format%tuple(head) + r'\\ \hline \hline' 
    for this_row in rows:
        head += '\n' + row_format%tuple(this_row) + r'\\'
    closing = "\\hline" + '\n'
    closing += r'\end{tabular} \end{center} \end{table}'
    return setup + head +closing
