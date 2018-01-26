from sage.all import cartesian_product, pi
from sage.misc.functional import numerical_approx as N
from sage.functions.other import sqrt
from sage.functions.trig import cos, sin

VERTEX = "v"
EDGE = "e"
FACE = "f"
LABELS = frozenset([VERTEX, EDGE, FACE])
DUAL = {VERTEX: FACE, EDGE: EDGE, FACE: VERTEX}
EDGE_COLORS = {VERTEX: "red", EDGE: "green", FACE: "blue"}
NONSIMPLE = {"loops": True, "multiedges": True,
             "immutable": True, "format": 'list_of_edges'}

T = sqrt(2 - sqrt(2)) / 2
S2 = N(sqrt(2))
S3 = N(sqrt(3)/2)
C8 = N(cos(pi/8) * T)
S8 = N(sin(pi/8) * T)
R = N(2*T*S3 + C8/T)

HEXAHEDRON_POS = {v: (N(r*cos(x)), N(r*sin(x))) for v, (r, x)
                  in enumerate((1+i//4, (2*i+1)*pi/4) for i in range(8))}

HORIZONTAL_LABELS = ['r', 'l']
VERTICAL_LABELS = ['u', 'd']

HEMIHEXAHEDRAL_FACES = [[1, 0, 3, 2], [2, 0, 1, 3], [3, 0, 2, 1]]

DODECAGON_FLAGS = ['au', 'Bu', 'Bd', 'cd', 'cu', 'Au',
                   'Ad', 'bd', 'bu', 'Cu', 'Cd', 'ad']
DODECAGON = {tuple(k): (N(cos(x) * T), N(sin(x) * T))
             for k, x in zip(DODECAGON_FLAGS,
                             [(2*i+1)*pi/12 for i in range(12)])}
DODECAGON_WRAP = {tuple(k): p for k, p
                  in zip(DODECAGON_FLAGS,
                         [(False, True), (False, True), (True, True),
                          (True, True), (True, True), (True, True),
                          (True, False), (True, False), (False, False),
                          (False, False), (False, False), (False, False)])}

TRIANGULAR_ARCS = 'abcBCA'
DODECAGON_SWAP = dict(zip(TRIANGULAR_ARCS, reversed(TRIANGULAR_ARCS)))

DODECAGON2_FLAGS = ['cl', 'cr', 'Ar', 'Al', 'bl', 'br',
                    'Cr', 'Cl', 'al', 'ar', 'Br', 'Bl']
DODECAGON2 = {tuple(k): (N(cos(x) * T), N(sin(x) * T))
              for k, x in zip(DODECAGON2_FLAGS,
                              [(2*i+1)*pi/12 for i in range(12)])}
DODECAGON2_WRAP = {tuple(k): p for k, p
                   in zip(DODECAGON2_FLAGS,
                          [(False, 1), (False, 2), (False, 2),
                           (True, 2), (True, 2), (True, 1),
                           (True, 1), (True, 0), (True, 0),
                           (False, 0), (False, 0), (False, 1)])}
DODECAGON2_SWAP = dict(zip('abcABC', 'acbACB'))

HORIZONTAL_SWAP = dict(zip(HORIZONTAL_LABELS, reversed(HORIZONTAL_LABELS)) +
                       zip(VERTICAL_LABELS, VERTICAL_LABELS))
HORIZONTAL_OFFSET = dict(zip(HORIZONTAL_LABELS, [1, -1]))
SQUARE_FLAGS = cartesian_product([HORIZONTAL_LABELS, VERTICAL_LABELS])
OCTAGON = {tuple(k): (N(cos(x) * T), N(sin(x) * T))
           for k, x in zip(['ru', 'ur', 'ul', 'lu', 'ld', 'dl', 'dr', 'rd'],
                           [(2*i+1)*pi/8 for i in range(8)])}
