from sage.all import cartesian_product, Integers
from sage.rings.integer import Integer
from .constants import VERTEX, EDGE, FACE, DODECAGON, OCTAGON
from .constants import TRIANGULAR_ARCS, TRIANGULAR_FACES
from .constants import HORIZONTAL_ARCS, VERTICAL_ARCS, SQUARE_FLAGS
from functions import squareEdgeFunction
from functions import squareFaceFunction
from functions import squarePosition
from functions import triangularEdgeFunction
from functions import triangularFaceFunction
from functions import triangularPosition
from functions import vertexFunction
from .tiling import Tiling

class Surface:
    def __init__(self):
        raise NotImplementedError

class Torus(Surface):
    @staticmethod
    def hexagonalTiling(n, m, k, center_faces = True):
        return Torus.triangularTiling(n, m, k,
                                      center_faces = center_faces).dual()

    @staticmethod
    def squareTiling(n, m, k):
        k = k % n
        vertices = cartesian_product([Integers(n), Integers(m)])
        edges = []
        pos = {}
        r, l = HORIZONTAL_ARCS
        u, d = VERTICAL_ARCS
        for v in vertices:
            v = tuple(v)
            i, j = v
            ib = i+k if j == -1 else i
            for hor in HORIZONTAL_ARCS:
                edges.append(((v, d, hor), ((ib, j+1), u, hor), VERTEX))
                edges.append(((v, hor, u), (v, hor, d), FACE))
            for ver in VERTICAL_ARCS:
                edges.append(((v, r, ver), ((i+1, j), l, ver), VERTEX))
                edges.append(((v, ver, l), (v, ver, r), FACE))
            for hor, ver in SQUARE_FLAGS:
                edges.append(((v, hor, ver), (v, ver, hor), EDGE))
            for (e, f), p in OCTAGON.items():
                pos[v, e, f] = [sum(p) for p in zip(squarePosition(v), p)]
        return Tiling(edges, pos = pos,
                      vertex_fun = vertexFunction,
                      edge_fun = squareEdgeFunction(k),
                      face_fun = squareFaceFunction(k))

    @staticmethod
    def triangularTiling(n, m, k, center_faces = False):
        k = k % n
        vertices = cartesian_product([Integers(n), Integers(m)])
        edges = []
        pos = {}
        wrap = (n, m, k) if center_faces else None
        for v in vertices:
            v = tuple(v)
            i, j = v
            ib = i+k if j == -1 else i
            ic = i-k if j == 0 else i
            edges.append(((v, 'a', 'u'), (v, 'B', 'u'), EDGE))
            edges.append(((v, 'a', 'd'), (v, 'C', 'd'), EDGE))
            edges.append(((v, 'b', 'u'), (v, 'C', 'u'), EDGE))
            edges.append(((v, 'b', 'd'), (v, 'A', 'd'), EDGE))
            edges.append(((v, 'c', 'u'), (v, 'A', 'u'), EDGE))
            edges.append(((v, 'c', 'd'), (v, 'B', 'd'), EDGE))
            for f in TRIANGULAR_FACES:
                edges.append(((v, 'a', f), ((i+1, j), 'A', f), VERTEX))
                edges.append(((v, 'b', f), ((ib, j+1), 'B', f), VERTEX))
                edges.append(((v, 'c', f), ((ic-1, j-1), 'C', f), VERTEX))
            for a in TRIANGULAR_ARCS:
                edges.append(((v, a, 'u'), (v, a, 'd'), FACE))
                for f in TRIANGULAR_FACES:
                    pos[v, a, f] = [sum(p) for p
                                    in zip(triangularPosition(v, a, f, wrap),
                                           DODECAGON[a, f])]
        return Tiling(edges, pos = pos,
                      vertex_fun = vertexFunction,
                      edge_fun = triangularEdgeFunction(k),
                      face_fun = triangularFaceFunction(k))
