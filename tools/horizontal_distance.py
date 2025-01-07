from starter2 import *

#stolen from here https://stackoverflow.com/questions/68695541/how-to-calculate-horizontal-distance-between-two-lines

def pairwise(iterable):
    from itertools import tee

    a, b = tee(iterable)
    _ = next(b, None)
    yield from zip(a, b)
class rake():
    def __init__(self,a,b):
        intersections = []

        tarr = np.zeros([len(a),len(b)])-2
        uarr = np.zeros([len(a),len(b)])-2
        pxarr = np.zeros(len(a))-2
        pyarr = np.zeros(len(a))-2
        tgood = np.zeros(len(a))-2
        for x1, y1 in enumerate(a):
            x2 = len(b)
            y2 = y1
            for x3, (y3, y4) in enumerate(pairwise(b)):
                x4 = x3 + 1

                try:
                    t = ((x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4)) / ((x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4))

                    u = ((x2 - x1) * (y1 - y3) - (y2 - y1) * (x1 - x3)) / ((x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4))
                except ZeroDivisionError:
                    continue

                tarr[x1,x3]=t
                uarr[x1,x3]=u
                if 0 <= t <= 1.0 and 0 <= u <= 1.0:
                    px, py = x1 + t * (x2 - x1), y1 + t * (y2 - y1)
                    intersections.append((x1, px, py))
                    pxarr[x1]=px
                    pyarr[x1]=py
                    tgood[x1]=t
                    break
        self.tgood=tgood
        self.y1 = copy.copy(a)
        self.y3 = copy.copy(b)
        self.uarr=uarr
        self.tarr=tarr
        self.pxarr=pxarr
        self.pyarr=pyarr



def ho2(xa=None,xb=None,ya=None,yb=None):
    y1 = nar(copy.copy(ya))
    y2 = nar(copy.copy(y1)) #horizontal lines
    y3 = nar(copy.copy(yb))
    y4 = copy.copy(yb) #will be trimmed to shift later
    if xa is None:
        xa = np.arange(y1.size)
    if xb is None:
        xb = np.arange(y3.size)
    intersections = []
    x1 = nar(copy.copy(xa))
    x2 = len(yb)+np.zeros_like(x1) #horizontal line to the end
    x3 = nar(copy.copy(xb))
    x4 = copy.copy(x3) #will get trimmed to shift 

    #trim
    y1=y1[:-1]
    y2=y2[:-1]
    x1=x1[:-1]
    x2=x2[:-1]
    y3=y3[:-1]
    y4=y4[1:] #here we shift.
    x3=x3[:-1]
    x4=x4[1:]

    ones_x = np.ones([y1.size])
    ones_y = np.ones([y3.size,1])


    x3.shape = x3.size,1
    x4.shape = x4.size,1
    y3.shape = y3.size,1
    y4.shape = y4.size,1

    y3 -= 1e-12
    y4 += 1e-8
    denom = ((x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4))
    t = np.zeros_like(x1-x3,dtype='float')-1
    u = np.zeros_like(x1-x3,dtype='float')-1
    ok = np.abs(denom)>0
    t[ok] =  ((x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4))[ok] / denom[ok]
    u[ok] = ((x2 - x1) * (y1 - y3) - (y2 - y1) * (x1 - x3))[ok] / denom[ok]
    ok2 = (0<=t)*(t<=1.0)*(0<=u)*(u<=1.0)

    #for each x1, we want the lowest x3.
    #AX1 and AX3 are all the segments that intersect.
    #unique gives the first instance of each in AX1
    #AX1,AX3 = np.where(ok2)
    #UX1,IX1 = np.unique(AX1,return_index=True)
    #UX3=AX3[IX1] #the lowest X3.
    AX1,AX3 = np.where(ok2)
    UX3,IX3 = np.unique(AX3, return_index=True)
    UX1 = AX1[IX3]
    pdb.set_trace()

    pxarr = x1 + t * (x2 - x1 )
    pyarr = y1 + t * (y2 - y1 )
    px = pxarr[(UX1,UX3)]
    py = pyarr[(UX1,UX3)]
    good_x1 = (x1*ones_y)[(UX1,UX3)]


    return np.stack([good_x1, px, py]).transpose()


def ho(a,b):
    intersections = []

    for x1, y1 in enumerate(a):
        x2 = len(b)
        y2 = y1
        for x3, (y3i, y4) in enumerate(pairwise(b)):
            x4 = x3 + 1
            y3 =y3i#-1e-3

            try:
                t = ((x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4)) / ((x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4))
                u = ((x2 - x1) * (y1 - y3) - (y2 - y1) * (x1 - x3)) / ((x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4))
            except ZeroDivisionError:
                continue

            if 0 <= t <= 1.0 and 0 <= u <= 1.0:
                px, py = x1 + t * (x2 - x1), y1 + t * (y2 - y1)
                intersections.append((x1, px, py))
                break
    return nar(intersections)

def ho3(a,b):
    intersections = []

    for x1, y1 in enumerate(a):
        x2 = len(b)
        y2 = y1
        for x3, (y3i, y4) in enumerate(pairwise(b)):
            x4 = x3 + 1
            y3 =y3i#-1e-3

            denom = ((x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4))
            if denom == 0:
                continue

            if np.abs(denom)>0:
                t = ((x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4)) / denom
                u = ((x2 - x1) * (y1 - y3) - (y2 - y1) * (x1 - x3)) / denom

            if 0 <= t <= 1.0 and 0 <= u <= 1.0:
                px, py = x1 + t * (x2 - x1), y1 + t * (y2 - y1)
                intersections.append((x1, px, py))
                break
        else:
            intersections.append((x1,x1,y1))
    return nar(intersections)

def try2(a,b,method=1,fname='hor_test_2'):
    if method==1:
        intersections = ho(a,b)
    if method==2:
        intersections = ho2(ya=a,yb=b)
    fig,axes=plt.subplots(1,2)
    ax0=axes[0];ax1=axes[1]
    ax0.plot(range(len(a)), a, color="blue")
    ax0.plot(range(len(b)), b, color="orange")
    dx=[]
    for ii, x, y in intersections:
        i=int(ii)
        xs = [i, x]
        dx.append(x-i)
        ys = [a[i], y]
        ax0.plot(xs, ys, "r--")
        #plt.plot(x, y, "r+")
    epb.equal_prob(nar(dx), 16, ax=ax1)
    fig.savefig('plots_to_sort/%s'%fname)
    return intersections
def try1(method=1, fname='hor_test_1'):
    a = [4, 1, 2, 7, 8, 8, 6, 11, 7, 10, 11, 15, 14, 14, 13, 17, 17, 21, 22, 20]
    b = [3, 0, 1, 6, 3, 6, 9, 11, 8, 8, 11, 15, 14, 15, 17, 14, 18, 17, 18, 20,21]
    if method==1:
        intersections = ho(a,b)
    if method==2:
        intersections = ho2(ya=a,yb=b)
    plt.clf()
    plt.plot(range(len(a)), a, color="blue")
    plt.plot(range(len(b)), b, color="orange")
    for i1, x, y in intersections:
        i = int(i1)
        xs = [i, x]
        ys = [a[i], y]
        plt.plot(xs, ys, "r--")
        plt.plot(x, y, "r+")
    plt.savefig('plots_to_sort/%s'%fname)
    return intersections
