from sage.all import cartesian_product, pi
from sage.misc.functional import numerical_approx as n
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

T = 3
S2 = n(sqrt(2))
S3 = n(sqrt(3)/2)
C8 = n(cos(pi/8)/T)
S8 = n(sin(pi/8)/T)

DODECAGON_FLAGS = ['au', 'Bu', 'Bd', 'cd', 'cu', 'Au',
                   'Ad', 'bd', 'bu', 'Cu', 'Cd', 'ad']
DODECAGON = {tuple(k): (n(cos(x)/T), n(sin(x)/T))
             for k, x in zip(DODECAGON_FLAGS,
                             [(2*i+1)*pi/12 for i in range(12)])}
DODECAGON_WRAP = {tuple(k): p for k, p
                  in zip(DODECAGON_FLAGS,
                         [(False, True), (False, True), (True, True),
                          (True, True), (True, True), (True, True),
                          (True, False), (True, False), (False, False),
                          (False, False), (False, False), (False, False)])}

TRIANGULAR_ARCS = 'abcABC'
TRIANGULAR_FACES = 'ud'

HORIZONTAL_ARCS = ['r', 'l']
VERTICAL_ARCS = ['u', 'd']
HORIZONTAL_SWAP = dict(zip(HORIZONTAL_ARCS, reversed(HORIZONTAL_ARCS)) +
                       zip(VERTICAL_ARCS, VERTICAL_ARCS))
HORIZONTAL_OFFSET = dict(zip(HORIZONTAL_ARCS, [1, -1]))
SQUARE_FLAGS = cartesian_product([HORIZONTAL_ARCS, VERTICAL_ARCS])
OCTAGON = {tuple(k): (n(cos(x)/T), n(sin(x)/T))
           for k, x in zip(['ru', 'ur', 'ul', 'lu', 'ld', 'dl', 'dr', 'rd'],
                           [(2*i+1)*pi/8 for i in range(8)])}
