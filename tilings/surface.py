from sage.all import cartesian_product, Integers, pi
from sage.functions.trig import cos, sin
from sage.graphs.graph_generators import graphs
from sage.misc.functional import numerical_approx as N
from sage.rings.integer import Integer
from .constants import R, S2, S3, VERTEX, EDGE, FACE
from .constants import HEXAHEDRON_POS, HEMIHEXAHEDRAL_FACES
from .constants import DODECAGON, DODECAGON2, OCTAGON
from .constants import HORIZONTAL_LABELS, VERTICAL_LABELS
from .constants import TRIANGULAR_ARCS, SQUARE_FLAGS
from .constants import HORIZONTAL_SWAP, HORIZONTAL_OFFSET
from functions import first, second
from functions import hosohedralFaceFunction
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

class Sphere(Surface):
    @staticmethod
    def dodecahedron():
        return Tiling(skeleton = graphs.DodecahedralGraph(),
                      faces = [[0, 1, 2, 3, 19], [0, 10, 11, 18, 19],
                               [0, 1, 8, 9, 10], [9, 10, 11, 12, 13],
                               [1, 2, 6, 7, 8], [11, 12, 16, 17, 18],
                               [2, 3, 4, 5, 6], [12, 13, 14, 15, 16],
                               [3, 4, 17, 18, 19], [5, 6, 7, 14, 15],
                               [4, 5, 15, 16, 17], [7, 8, 9, 13, 14]])

    greatStellatedDodecahedron = dodecahedron

    @staticmethod
    def greatDodecahedron():
        I = graphs.IcosahedralGraph()
        I._pos[3] = (15*S3, -7.5)
        I._pos[9] = (0, 15)
        I._pos[10] = (-15*S3, -7.5)
        return Tiling(skeleton = I,
                      faces = [[0, 1, 2, 9, 7], [0, 8, 9, 10, 11],
                               [0, 1, 6, 4, 11], [0, 5, 6, 2, 8],
                               [0, 5, 4, 10, 7], [3, 6, 5, 11, 10],
                               [1, 2, 3, 4, 5], [1, 5, 11, 7, 8],
                               [2, 3, 10, 7, 8], [2, 6, 4, 10, 9],
                               [1, 6, 3, 9, 8], [3, 4, 11, 7, 9]])

    @staticmethod
    def hexahedron():
        H = graphs.HexahedralGraph()
        H._pos = HEXAHEDRON_POS
        return Tiling(skeleton = H,
                      faces = [[0, 1, 2, 3], [0, 1, 5, 4], [0, 3, 7, 4],
                               [1, 2, 6, 5], [2, 3, 7, 6], [4, 5, 6, 7]])

    @staticmethod
    def hosohedron(n, center_faces = False):
        edges = []
        pos = {}
        pin = pi/n if center_faces else pi/(2*n-1)
        for e in Integers(n):
            for f in HORIZONTAL_LABELS:
                edges.append((('u', e, f), ('d', e, f), VERTEX))
            for v in VERTICAL_LABELS:
                edges.append(((v, e, 'r'), (v, e+1, 'l'), EDGE))
                edges.append(((v, e, 'r'), (v, e, 'l'), FACE))
            i = 2*Integer(e)
            if center_faces:
                pos['u', e, 'l'] = (N(n*cos(i*pin)), N(n*sin(i*pin)))
                pos['u', e, 'r'] = (N(n*cos((i+1)*pin)), N(n*sin((i+1)*pin)))
                pos['d', e, 'l'] = (N((n-1)*cos(i*pin)), N((n-1)*sin(i*pin)))
                pos['d', e, 'r'] = (N((n-1)*cos((i+1)*pin)),
                                    N((n-1)*sin((i+1)*pin)))
            else:
                pos['u', e, 'l'] = (N(-n*cos(i*pin)), N(n*(1-sin(i*pin))))
                pos['u', e, 'r'] = (N(-n*cos((i+1)*pin)),
                                    N(n*(1-sin((i+1)*pin))))
                pos['d', e, 'l'] = (N(-n*cos(i*pin)), N(n*(sin(i*pin)-1)))
                pos['d', e, 'r'] = (N(-n*cos((i+1)*pin)),
                                    N(n*(sin((i+1)*pin)-1)))
        return Tiling(edges, pos = pos,
                      vertex_fun = first,
                      edge_fun = second,
                      face_fun = hosohedralFaceFunction)

    @staticmethod
    def icosahedron():
        I = graphs.IcosahedralGraph()
        I._pos[3] = (15*S3, -7.5)
        I._pos[9] = (0, 15)
        I._pos[10] = (-15*S3, -7.5)
        return Tiling(skeleton = I,
                      faces = [[0, 1, 5], [0, 1, 8], [0, 5, 11], [0, 7, 11],
                               [0, 7, 8], [1, 2, 6], [1, 2, 8], [4, 10, 11],
                               [1, 5, 6], [2, 3, 6], [2, 3, 9], [7, 10, 11],
                               [2, 8, 9], [3, 4, 6], [3, 4, 10], [3, 9, 10],
                               [4, 5, 6], [4, 5, 11], [7, 8, 9], [7, 9, 10]])

    greatIcosahedron = icosahedron

    @staticmethod
    def octahedron():
        return Tiling(skeleton = graphs.OctahedralGraph(),
                      faces = [[0, 1, 2], [0, 1, 3], [0, 2, 4], [0, 3, 4],
                               [1, 2, 5], [1, 3, 5], [2, 4, 5], [3, 4, 5]])

    @staticmethod
    def polygon(n, center_faces = True):
        return Sphere.hosohedron(n, center_faces = center_faces).dual()

    @staticmethod
    def smallStellatedDodecahedron():
        return Sphere.greatDodecahedron().dual()

    @staticmethod
    def tetrahedron():
        return Tiling(skeleton = graphs.TetrahedralGraph(),
                      faces = [[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]])

class ProjectivePlane(Surface):
    @staticmethod
    def hemidodecahedron():
        return Tiling(skeleton = graphs.PetersenGraph(),
                      faces = [[0, 1, 2, 3, 4], [0, 1, 6, 8, 5],
                               [0, 4, 9, 7, 5], [1, 2, 7, 9, 6],
                               [2, 3, 8, 5, 7], [3, 4, 9, 6, 8]])

    @staticmethod
    def hemihexahedron():
        H = Tiling(skeleton = graphs.TetrahedralGraph(),
                   faces = HEMIHEXAHEDRAL_FACES)
        pos = {}
        for f, c in enumerate(HEMIHEXAHEDRAL_FACES):
            phi = 2*f*pi/3
            cx = N(R * cos(phi))
            cy = N(R * sin(phi))
            for i in range(4):
                pos[c[i-1], c[i], f] = (cx + cos(phi + (1+4*i)*pi/8),
                                        cy + sin(phi + (1+4*i)*pi/8))
                pos[c[i], c[i-1], f] = (cx + cos(phi + (3+4*i)*pi/8),
                                        cy + sin(phi + (3+4*i)*pi/8))
        H.set_pos(pos)
        return H

    @staticmethod
    def hemihosohedron(n):
        edges = []
        pos = {}
        phi = pi/2 if n == 1 else pi
        for e in Integers(n):
            i = 4*Integer(e)
            phl, phr = (i-1)*phi/(2*n), (i+1)*phi/(2*n)
            ru = n + (1 + sin(phr/2))/2
            lu = n + (1 + sin(phl/2))/2
            rd = n + (1 - sin(phr/2))/2
            ld = n + (1 - sin(phl/2))/2
            varcs, lu, ld = (reversed(VERTICAL_LABELS), ld, lu) if e == 0 \
                       else (VERTICAL_LABELS, lu, ld)
            for h in HORIZONTAL_LABELS:
                edges.append((('u', e, h), ('d', e, h), VERTEX))
            for v, ver in zip(VERTICAL_LABELS, varcs):
                edges.append(((v, e, 'r'), (v, e+1, 'l'), EDGE))
                edges.append(((v, e, 'r'), (ver, e, 'l'), FACE))
            pos['u', e, 'r'] = (N(ru*cos(phr)), N(ru*sin(phr)))
            pos['u', e, 'l'] = (N(ru*cos(phl)), N(ru*sin(phl)))
            pos['d', e, 'r'] = (N(rd*cos(phr)), N(rd*sin(phr)))
            pos['d', e, 'l'] = (N(rd*cos(phl)), N(rd*sin(phl)))
        return Tiling(edges, pos = pos,
                      vertex_fun = lambda x: (),
                      edge_fun = second,
                      face_fun = hosohedralFaceFunction)

    @staticmethod
    def hemiicosahedron():
        return ProjectivePlane.hemidodecahedron().dual()

    @staticmethod
    def hemioctahedron():
        return ProjectivePlane.hemihexahedron().dual()

    @staticmethod
    def hemipolygon(n):
        return ProjectivePlane.hemihosohedron(n).dual()

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
