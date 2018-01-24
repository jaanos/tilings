from sage.all import cartesian_product, Integers
from sage.rings.integer import Integer
from .constants import S2, S3, VERTEX, EDGE, FACE
from .constants import DODECAGON, DODECAGON2, OCTAGON
from .constants import HORIZONTAL_LABELS, VERTICAL_LABELS
from .constants import TRIANGULAR_ARCS, SQUARE_FLAGS
from .constants import HORIZONTAL_SWAP, HORIZONTAL_OFFSET
from functions import kleinBottleSquareEdgeFunction1
from functions import kleinBottleSquareFaceFunction1
from functions import kleinBottleSquarePosition1
from functions import kleinBottleSquareWrap1
from functions import kleinBottleSquareEdgeFunction2
from functions import kleinBottleSquareFaceFunction2
from functions import kleinBottleSquareWrap2
from functions import kleinBottleTriangularEdgeFunction1
from functions import kleinBottleTriangularFaceFunction1
from functions import kleinBottleTriangularPosition1
from functions import kleinBottleTriangularWrap1
from functions import kleinBottleTriangularEdgeFunction2
from functions import kleinBottleTriangularFaceFunction2
from functions import kleinBottleTriangularPosition2
from functions import kleinBottleTriangularWrap2
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
        r, l = HORIZONTAL_LABELS
        u, d = VERTICAL_LABELS
        for v in vertices:
            v = tuple(v)
            i, j = v
            ib = i+k if j == -1 else i
            for hor in HORIZONTAL_LABELS:
                edges.append(((v, d, hor), ((ib, j+1), u, hor), VERTEX))
                edges.append(((v, hor, u), (v, hor, d), FACE))
            for ver in VERTICAL_LABELS:
                edges.append(((v, r, ver), ((i+1, j), l, ver), VERTEX))
                edges.append(((v, ver, l), (v, ver, r), FACE))
            for hor, ver in SQUARE_FLAGS:
                edges.append(((v, hor, ver), (v, ver, hor), EDGE))
            for (e, f), q in OCTAGON.items():
                pos[v, e, f] = [sum(p) for p in zip(squarePosition(v), q)]
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
            for f in VERTICAL_LABELS:
                edges.append(((v, 'a', f), ((i+1, j), 'A', f), VERTEX))
                edges.append(((v, 'b', f), ((ib, j+1), 'B', f), VERTEX))
                edges.append(((v, 'c', f), ((ic-1, j-1), 'C', f), VERTEX))
            for a in TRIANGULAR_ARCS:
                edges.append(((v, a, 'u'), (v, a, 'd'), FACE))
                for f in VERTICAL_LABELS:
                    pos[v, a, f] = [sum(p) for p
                                    in zip(triangularPosition(v, a, f, wrap),
                                           DODECAGON[a, f])]
        return Tiling(edges, pos = pos,
                      vertex_fun = vertexFunction,
                      edge_fun = triangularEdgeFunction(k),
                      face_fun = triangularFaceFunction(k))

class KleinBottle(Surface):
    @staticmethod
    def hexagonalTiling1(n, m, center_faces = True):
        return KleinBottle.triangularTiling1(n, m,
                                        center_faces = center_faces).dual()

    @staticmethod
    def hexagonalTiling2(n, m, center_faces = True):
        return KleinBottle.triangularTiling2(n, m,
                                        center_faces = center_faces).dual()

    @staticmethod
    def squareTiling1(n, m, center_faces = False):
        vertices = cartesian_product([Integers(n), Integers(m)])
        edges = []
        pos = {}
        wrap = (n, m) if center_faces else None
        r, l = HORIZONTAL_LABELS
        u, d = VERTICAL_LABELS
        for v in vertices:
            v = tuple(v)
            i, j = v
            ib, harcs = (-i, reversed(HORIZONTAL_LABELS)) if j == -1 \
                        else (i, HORIZONTAL_LABELS)
            for hor, horb in zip(HORIZONTAL_LABELS, harcs):
                edges.append(((v, d, hor), ((ib, j+1), u, horb), VERTEX))
                edges.append(((v, hor, u), (v, hor, d), FACE))
            for ver in VERTICAL_LABELS:
                edges.append(((v, r, ver), ((i+1, j), l, ver), VERTEX))
                edges.append(((v, ver, l), (v, ver, r), FACE))
            for hor, ver in SQUARE_FLAGS:
                edges.append(((v, hor, ver), (v, ver, hor), EDGE))
            for (e, f), p in OCTAGON.items():
                pos[v, e, f] = [sum(p) for p
                                in zip(kleinBottleSquarePosition1(v, e, f,
                                                                  wrap),
                                       kleinBottleSquareWrap1(p, v, e, f,
                                                              wrap))]
        return Tiling(edges, pos = pos,
                      vertex_fun = vertexFunction,
                      edge_fun = kleinBottleSquareEdgeFunction1,
                      face_fun = kleinBottleSquareFaceFunction1)

    @staticmethod
    def squareTiling1Dual(n, m, center_faces = True):
        return KleinBottle.squareTiling1(n, m,
                                         center_faces = center_faces).dual()

    @staticmethod
    def squareTiling2(n, m):
        vertices = [(i, j) for i, j in cartesian_product([Integers(2*n),
                                                          Integers(m)])
                    if (Integer(i) + Integer(j)) % 2 == 0]
        edges = []
        pos = {}
        r, l = HORIZONTAL_LABELS
        u, d = VERTICAL_LABELS
        for v in vertices:
            i, j = v
            for hor in HORIZONTAL_LABELS:
                ib, horb = kleinBottleSquareWrap2(i, j+1, hor, m)
                horb = HORIZONTAL_SWAP[horb]
                ib -= HORIZONTAL_OFFSET[horb]
                edges.append(((v, d, hor), ((ib, j+1), horb, u), VERTEX))
                edges.append(((v, hor, d), ((ib, j+1), u, horb), VERTEX))
                edges.append(((v, hor, u), (v, hor, d), EDGE))
            for ver in VERTICAL_LABELS:
                edges.append(((v, ver, l), (v, ver, r), EDGE))
            for hor, ver in SQUARE_FLAGS:
                edges.append(((v, hor, ver), (v, ver, hor), FACE))
            for (e, f), q in OCTAGON.items():
                pos[v, e, f] = [x + y * S2 for x, y
                                in zip(squarePosition(v), q)]
        return Tiling(edges, pos = pos,
                      vertex_fun = vertexFunction,
                      edge_fun = kleinBottleSquareEdgeFunction2(m),
                      face_fun = kleinBottleSquareFaceFunction2(m))

    @staticmethod
    def triangularTiling1(n, m, center_faces = False):
        vertices = cartesian_product([Integers(n), Integers(m)])
        edges = []
        pos = {}
        wrap = (n, m) if center_faces else None
        for v in vertices:
            v = tuple(v)
            i, j = v
            ib, ic, bb, cc = (-i, -i-1, 'c', 'B') if j == -1 else \
                             (i, i+1, 'B', 'c')
            edges.append(((v, 'a', 'u'), (v, 'B', 'u'), EDGE))
            edges.append(((v, 'a', 'd'), (v, 'C', 'd'), EDGE))
            edges.append(((v, 'b', 'u'), (v, 'C', 'u'), EDGE))
            edges.append(((v, 'b', 'd'), (v, 'A', 'd'), EDGE))
            edges.append(((v, 'c', 'u'), (v, 'A', 'u'), EDGE))
            edges.append(((v, 'c', 'd'), (v, 'B', 'd'), EDGE))
            for f in VERTICAL_LABELS:
                edges.append(((v, 'a', f), ((i+1, j), 'A', f), VERTEX))
                edges.append(((v, 'b', f), ((ib, j+1), bb, f), VERTEX))
                edges.append(((v, 'C', f), ((ic, j+1), cc, f), VERTEX))
            for a in TRIANGULAR_ARCS:
                edges.append(((v, a, 'u'), (v, a, 'd'), FACE))
            for (e, f), p in DODECAGON.items():
                pos[v, e, f] = [sum(p) for p
                                in zip(kleinBottleTriangularPosition1(v, e, f,
                                                                      wrap),
                                       kleinBottleTriangularWrap1(p, v, e, f,
                                                                  wrap))]
        return Tiling(edges, pos = pos,
                      vertex_fun = vertexFunction,
                      edge_fun = kleinBottleTriangularEdgeFunction1,
                      face_fun = kleinBottleTriangularFaceFunction1)

    @staticmethod
    def triangularTiling2(n, m, center_faces = False):
        vertices = [(i, j) for i, j in cartesian_product([Integers(2*n),
                                                          Integers(m)])
                    if (Integer(i) + Integer(j)) % 2 == 0]
        edges = []
        pos = {}
        wrap = (n, m) if center_faces else None
        h = m % 2
        for v in vertices:
            v = tuple(v)
            i, j = v
            ia, hfa = (-i-h, reversed(HORIZONTAL_LABELS)) if j in [-1, -2] \
                 else (i, HORIZONTAL_LABELS)
            ib, ic, bb, cc, hfb = (-i-1-h, -i+1-h, 'c', 'b',
                                   reversed(HORIZONTAL_LABELS)) \
                                  if j == -1 else \
                                  (i+1, i-1, 'b', 'c', HORIZONTAL_LABELS)
            edges.append(((v, 'a', 'r'), (v, 'B', 'r'), EDGE))
            edges.append(((v, 'a', 'l'), (v, 'C', 'l'), EDGE))
            edges.append(((v, 'b', 'r'), (v, 'C', 'r'), EDGE))
            edges.append(((v, 'b', 'l'), (v, 'A', 'l'), EDGE))
            edges.append(((v, 'c', 'r'), (v, 'A', 'r'), EDGE))
            edges.append(((v, 'c', 'l'), (v, 'B', 'l'), EDGE))
            for hor, hora, horb in zip(HORIZONTAL_LABELS, hfa, hfb):
                edges.append(((v, 'a', hor), ((ia, j+2), 'A', hora), VERTEX))
                edges.append(((v, 'B', hor), ((ib, j+1), bb, horb), VERTEX))
                edges.append(((v, 'C', hor), ((ic, j+1), cc, horb), VERTEX))
            for a in TRIANGULAR_ARCS:
                edges.append(((v, a, 'l'), (v, a, 'r'), FACE))
            for (e, f), p in DODECAGON2.items():
                pos[v, e, f] = [x + y / S3 for x, y
                                in zip(kleinBottleTriangularPosition2(v, e, f,
                                                                      wrap),
                                       kleinBottleTriangularWrap2(p, v, e, f,
                                                                  wrap))]
        return Tiling(edges, pos = pos,
                      vertex_fun = vertexFunction,
                      edge_fun = kleinBottleTriangularEdgeFunction2(m),
                      face_fun = kleinBottleTriangularFaceFunction2(m))
