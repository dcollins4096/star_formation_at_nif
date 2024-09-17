
from starter2 import *
import yt
class gitter():
    def __init__(self,template):
        self.template=template
    def __call__(self,frame,field):
        self.ds = yt.load( self.template%(frame,frame))
        ad = self.ds.all_data()
        x = ad['x']
        order = np.argsort(x)
        field = ad[field][order]
        return field
class gitter2():
    def __init__(self,template):
        self.template=template
    def __call__(self,frame,field):
        self.ds = yt.load( self.template%(frame,frame))
        ray = self.ds.ortho_ray(0,[0,0])
        field = ray[field]
        return field

