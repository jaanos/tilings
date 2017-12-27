from sage.rings.integer import Integer
from .constants import S3

def meanpos(G, l):
    return tuple(sum(p) / float(len(p))
                 for p in zip(*(G._pos[x] for x in l)))

triangularPosition = lambda (i, j): (Integer(i) - Integer(j)/2,
                                     -Integer(j) * S3)

triangularVertexFunction = lambda (v, e, f): v

def triangularEdgeFunction(k):
    def edgeFun((v, e, f)):
        i, j = v
        if e == 'A':
            return ((i-1, j), 'a')
        elif e == 'B':
            return ((i-k if j == 0 else i, j-1), 'b')
        elif e == 'C':
            return ((i+k+1 if j == -1 else i+1, j+1), 'c')
        else:
            return (v, e)
    return edgeFun

def triangularFaceFunction(k):
    def faceFun((v, e, f)):
        i, j = v
        u = v
        if f == 'u':
            if e in 'aB':
                u = (i-k if j == 0 else i, j-1)
            elif e in 'Ac':
                u = (i-k-1 if j == 0 else i-1, j-1)
        else:
            if e in 'Ab':
                u = (i+k if j == -1 else i, j+1)
            elif e in 'aC':
                u = (i+k+1 if j == -1 else i+1, j+1)
        return (u, f)
    return faceFun