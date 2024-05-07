from starter2 import *


def frmt(value, f="%0.2f", unit=None):
    un = value.units
    if unit is not None:
        un = unit
    return (f+" %s")%(value.v, un)

def table(devices,fname):

    head = ['shot',r'$I_{01}$','$I_{02}$','Z',r'$\Delta x$',r'$\Delta \rho$', r'$v_s$', r'$R$', r'$c_s$', r'$\bar{A}$', r'$\sigma_B$']
    head += [r'$\sigma_v$', r'$\mathcal{M}_s$']
    rows=[]
    kms=r' $\rm{km}\ \rm{s}^{-1}$'
    mgpercc=r' $\rm{mg}\ \rm{cm}^{-3}$'
    gpercc=r' $\rm{g}\ \rm{cm}^{-3}$'
    for name in devices:
        dev = devices[name]
        stuff = [name,
                 "%0.1f"%dev.shot1.I0,
                 "%0.1f"%dev.shot2.I0,
                 "%0.1f"%dev.shot1.zero_value,
                 frmt(dev.shift_x.in_units('um'), unit= r'$\mu m$'),
                 frmt(dev.shift_rho.in_units('mg/cm**3'),unit= mgpercc),
                 "%0.1f"%dev.vel.in_units('km/s') + kms,
                 "%0.1f"%dev.R,
                 "%0.1f"%dev.cs.in_units('km/s') + kms,
                 "%0.1f"%dev.atwood_number,
                 "%0.1f"%dev.sigma_B + gpercc,
                 "%0.1f"%dev.sigma_v.in_units('km/s') + kms,
                 "%0.1f"%(dev.sigma_v/dev.cs)
                ]
        this_row= "%s &"*len(stuff)%tuple(stuff)
        this_row=this_row[:-1] + "\\\\\n"
        rows.append(this_row)
    aligner = '{'+'r'*len(head)+'}'
    headline = '%s &'*len(head)%tuple(head) 
    headline = headline[:-1] #too many &
    headline += '\\\\\n' #yes a lot of backslash.
    hline=r'\hline'+'\n'
    fptr=open('%s/%s'%(plot_dir,fname),'w')
    fptr.write(r'\begin{table*}\begin{center}\input{table2_caption}\begin{tabular}'+aligner+'\n')
    fptr.write(hline)
    fptr.write(headline)
    fptr.write(hline)
    for row in rows:
        fptr.write(row)
    fptr.write(hline)
    fptr.write('\end{tabular}\end{center}\end{table*}\n')
    fptr.close()
        

