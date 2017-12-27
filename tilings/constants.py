from sage.all import pi
from sage.misc.functional import numerical_approx as n
from sage.functions.other import sqrt
from sage.functions.trig import cos, sin

VERTEX = "v"
EDGE = "e"
FACE = "f"
LABELS = frozenset([VERTEX, EDGE, FACE])
DUAL = {VERTEX: FACE, EDGE: EDGE, FACE: VERTEX}
EDGE_COLORS = {VERTEX: "red", EDGE: "green", FACE: "blue"}
NONSIMPLE = {"loops": True, "multiedges": True, "immutable": True}

T = 3
S3 = n(sqrt(3)/2)
C8 = n(cos(pi/8)/T)
S8 = n(sin(pi/8)/T)
DODECAGON = {tuple(k): (n(cos(x)/T), n(sin(x)/T))
             for k, x in zip(['au', 'Bu', 'Bd', 'cd', 'cu', 'Au',
                              'Ad', 'bd', 'bu', 'Cu', 'Cd', 'ad'],
                             [(2*i+1)*pi/12 for i in range(12)])}

TRIANGULAR_ARCS = 'abcABC'
TRIANGULAR_FACES = 'ud'