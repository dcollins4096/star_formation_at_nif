#
# Need a volume average tool.
#
import numpy as na
nar = na.array
def volavg(array, rank=None, refine_by=None,debug=0):
    """volume averages *array* by *refine_by* in each direction.
    *rank* is an explicit argument to keep the user from being confused when he/she 
    accidentally passes in a flattened array."""
    #rank = len(array.shape)
    if refine_by == 1:
        return array
    if rank == None:
        print("Rank must be done by hand, to prevent user stupidity.")
        print("I'll now raise a stupid error")
        raise
    if (na.array(array.shape) % refine_by != 0 ).any():
        print("Cannot downsample array of size", array.shape, "by ", refine_by)
        if debug:
            print(na.array(array.shape) % refine_by)
        return array
    last = array 
    all = slice(None,None,None)
    sub = [slice(i,array.shape[0]-(refine_by-i-1),refine_by) for i in range(refine_by)]
    sub[-1] = slice(refine_by-1,None,refine_by) #not elegant, must be a better way.
    all_set = [all]*rank
    r_inv = 1./refine_by
    for d in range(rank):
        next_dims = nar(last.shape)
        next_dims[d] /= refine_by
        next = na.zeros(next_dims)
        this_set = all_set[:]
        if debug:
            print("dims: last",last.shape, "next",next.shape)
        for i in range(refine_by):
            this_set[d] = sub[i]
            next += last[tuple(this_set)] * r_inv
        last = next
    return next  
if 0:
    test = na.arange(8)
    test.shape = (2,2,2)
    x = volavg(test,debug=1)
    

