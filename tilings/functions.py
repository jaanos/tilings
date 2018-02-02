from sage.graphs.graph import Graph
from sage.rings.integer import Integer
from sage.sets.set import Set
from .constants import VERTEX, EDGE, FACE, CORNER, S3, S8, C8
from .constants import DODECAGON_WRAP, DODECAGON_SWAP
from .constants import DODECAGON2_WRAP, DODECAGON2_SWAP
from .constants import HORIZONTAL_SWAP, HORIZONTAL_OFFSET
from .constants import HORIZONTAL_LABELS, VERTICAL_LABELS

makeEdge = lambda e: e*2 if len(e) == 1 else e

first = lambda x: x[0]

second = lambda x: x[1]

def simplestGraph(edges):
    G = Graph(edges, loops = True, multiedges = True, immutable = True,
              format = 'list_of_edges')
    loops = G.has_loops()
    multiedges = G.has_multiple_edges()
    if not (loops and multiedges):
        G = Graph(G, immutable = True, loops = loops, multiedges = multiedges)
    return G

def meanpos(G, l):
    return tuple(sum(p) / float(len(p))
                 for p in zip(*(G._pos[x] for x in l)))

flagPosition = lambda pu, pv, pw: tuple((1-C8-S8) * a +
                                        C8 * b + S8 * c
                                        for a, b, c in zip(pu, pv, pw))

def kleinBottleSquarePosition1(v, a, f, wrap = None):
    i, j = v
    if wrap is not None:
        n, m = wrap
        if 'u' in [a, f] and j == 0:
            i = -i
            j = m
            a = HORIZONTAL_SWAP[a]
            f = HORIZONTAL_SWAP[f]
        if 'l' in [a, f] and i == 0:
            i = n
    return (Integer(i), -Integer(j))

def kleinBottleSquareWrap1(p, v, a, f, wrap = None):
    i, j = v
    x, y = p
    return (-x, y) if wrap is not None and 'u' in [a, f] and j == 0 else p

kleinBottleSquareWrap2 = lambda i, j, hor, m: \
    (-i if m % 2 == 0 else -i-1, HORIZONTAL_SWAP[hor]) if j == 0 else (i, hor)

def squarePosition((i, j)):
    return (Integer(i), -Integer(j))

def kleinBottleTriangularPosition1(v, a, f, wrap = None):
    i, j = v
    if wrap is not None:
        n, m = wrap
        h, d = DODECAGON_WRAP[a, f]
        if d and j == 0:
            j = m
            i = -i
            a = DODECAGON_SWAP[a]
            h, d = DODECAGON_WRAP[a, f]
        if h and i == 0:
            i = n
    i = Integer(i)
    j = Integer(j)
    return (i - j/2, -j * S3)

def kleinBottleTriangularWrap1(p, v, a, f, wrap = None):
    i, j = v
    x, y = p
    return (-x, y) if wrap is not None and j == 0 and \
                      (a in 'Bc' or (f == 'u' and a in 'aA')) else p

def kleinBottleTriangularPosition2(v, a, f, wrap = None):
    i, j = v
    j = Integer(j)
    if wrap is not None:
        n, m = wrap
        k = m % 2
        h, d = DODECAGON2_WRAP[a, f]
        if d > j:
            j += m
            i = -i - k
            h, d = DODECAGON2_WRAP[DODECAGON2_SWAP[a], HORIZONTAL_SWAP[f]]
        if h and i == 0:
            i = 2*n
    i = Integer(i)
    return (i, -j / (2*S3))

def kleinBottleTriangularWrap2(p, v, a, f, wrap = None):
    i, j = v
    x, y = p
    return (-x, y) if wrap is not None and \
        ((j == 0 and (a in 'Abc' or ((a, f) in [('B', 'l'), ('C', 'r')]))) or
         (j == 1 and (a == 'A' or ((a, f) in [('b', 'l'), ('c', 'r')])))) \
        else p

def triangularPosition(v, a, f, wrap = None):
    i, j = v
    if wrap is not None:
        n, m, k = wrap
        h, d = DODECAGON_WRAP[a, f]
        if d and j == 0:
            j = m
            i -= k
        if h and i == 0:
            i = n
    i = Integer(i)
    j = Integer(j)
    return (i - j/2, -j * S3)

hosohedralFaceFunction = lambda x: x[-2] if x[-1] == 'r' else x[-2]-1

def kleinBottleSquareEdgeFunction1((v, e, f)):
    i, j = v
    if e == 'l':
        return ((i-1, j), 'r')
    elif e == 'u':
        return ((-i if j == 0 else i, j-1), 'd')
    else:
        return (v, e)

def kleinBottleSquareFaceFunction1((v, e, f)):
    i, j = v
    lb = 'l'
    if 'u' in [e, f]:
        if j == 0:
            i = -i
            lb = 'r'
        j -= 1
    if lb in [e, f]:
        i -= 1
    return (i, j)

def kleinBottleSquareEdgeFunction2(m):
    def edgeFun((v, e, f)):
        i, j = v
        hor = next(x for x in [e, f] if x in HORIZONTAL_LABELS)
        ver = next(x for x in [e, f] if x in VERTICAL_LABELS)
        if ver == 'u':
            i, hor = kleinBottleSquareWrap2(i, j, hor, m)
            i, j = (i + HORIZONTAL_OFFSET[hor], j-1)
            hor = HORIZONTAL_SWAP[hor]
        return ((i, j), hor)
    return edgeFun

def kleinBottleSquareFaceFunction2(m):
    def faceFun((v, e, f)):
        i, j = v
        if e == 'u':
            i, f = kleinBottleSquareWrap2(i, j, f, m)
            i, j = (i + HORIZONTAL_OFFSET[f], j-1)
            e, f = HORIZONTAL_SWAP[f], e
        if e in HORIZONTAL_LABELS:
            i, e = kleinBottleSquareWrap2(i, j, e, m)
            i, j = (i + HORIZONTAL_OFFSET[e], j-1)
        return (i, j)
    return faceFun

def squareEdgeFunction(k):
    def edgeFun((v, e, f)):
        i, j = v
        if e == 'l':
            return ((i-1, j), 'r')
        elif e == 'u':
            return ((i-k if j == 0 else i, j-1), 'd')
        else:
            return (v, e)
    return edgeFun

def squareFaceFunction(k):
    def faceFun((v, e, f)):
        i, j = v
        if 'l' in [e, f]:
            i -= 1
        if 'u' in [e, f]:
            if j == 0:
                i -= k
            j -= 1
        return (i, j)
    return faceFun

def kleinBottleTriangularEdgeFunction1((v, e, f)):
    i, j = v
    ib, bb, cc = (-i, 'c', 'B') if j == 0 else (i, 'B', 'c')
    if e == 'A':
        return ((i-1, j), 'a')
    elif e == bb:
        return ((ib, j-1), 'b')
    elif e == cc:
        return ((ib-1, j-1), 'C')
    else:
        return (v, e)

def kleinBottleTriangularFaceFunction1((v, e, f)):
    i, j = v
    if f == 'u':
        if e in 'aB':
            v = (-i-1 if j == 0 else i, j-1)
        elif e in 'Ac':
            v = (-i if j == 0 else i-1, j-1)
    else:
        if e in 'Ab':
            v = (-i if j == -1 else i, j+1)
        elif e in 'aC':
            v = (-i-1 if j == -1 else i+1, j+1)
    return (v, f)

def kleinBottleTriangularEdgeFunction2(m):
    h = m % 2
    def edgeFun((v, e, f)):
        i, j = v
        ib, bb, cc = (-i-h, 'c', 'b') if j == 0 else (i, 'b', 'c')
        if e == 'A':
            return ((-i-h if j in [0, 1] else i, j-2), 'a')
        elif e == bb:
            return ((ib-1, j-1), 'B')
        elif e == cc:
            return ((ib+1, j-1), 'C')
        else:
            return (v, e)
    return edgeFun

def kleinBottleTriangularFaceFunction2(m):
    h = m % 2
    def faceFun((v, e, f)):
        i, j = v
        if f == 'r':
            if e in 'aB':
                v = (-i-h-1 if j == -1 else i+1, j+1)
            elif e in 'Ac':
                v = (-i-h+1 if j == 0 else i+1, j-1)
        else:
            if e in 'Ab':
                v = (-i-h-1 if j == 0 else i-1, j-1)
            elif e in 'aC':
                v = (-i-h+1 if j == -1 else i-1, j+1)
        return (v, f)
    return faceFun

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
        if f == 'u':
            if e in 'aB':
                v = (i-k if j == 0 else i, j-1)
            elif e in 'Ac':
                v = (i-k-1 if j == 0 else i-1, j-1)
        else:
            if e in 'Ab':
                v = (i+k if j == -1 else i, j+1)
            elif e in 'aC':
                v = (i+k+1 if j == -1 else i+1, j+1)
        return (v, f)
    return faceFun

truncationVertexFunction = lambda s: Set(p for p, l in s)

vertexFunction = lambda (v, e, f): v
