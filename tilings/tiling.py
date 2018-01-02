from sage.graphs.graph import Graph
from .constants import VERTEX, EDGE, FACE, LABELS, DUAL, EDGE_COLORS
from .constants import NONSIMPLE, C8, S8
from .functions import meanpos

class Tiling(Graph):
    def __init__(self, *largs, **kargs):
        self._skeleton = None
        self._muscles = None
        self._vertex_fun = None
        self._edge_fun = None
        self._face_fun = None
        self._dual = kargs.pop("dual", None)
        vertex_fun = kargs.pop("vertex_fun", None)
        edge_fun = kargs.pop("edge_fun", None)
        face_fun = kargs.pop("face_fun", None)
        if self._dual is not None:
            self._skeleton = self._dual._muscles
            self._muscles = self._dual._skeleton
            self._vertex_fun = self._dual._face_fun
            self._edge_fun = self._dual._edge_fun
            self._face_fun = self._dual._vertex_fun
        faces = None
        if "skeleton" in kargs or "faces" in kargs:
            assert "skeleton" in kargs and "faces" in kargs, \
                "both skeleton and faces should be given, or none"
            assert len(largs) == 0 and kargs.get("data") is None, \
                "graph data should not be given along skeleton and faces"
            self._skeleton = kargs.pop("skeleton")
            faces = kargs.pop("faces")
            if not isinstance(faces, dict):
                faces = dict(enumerate(faces))
            self._skeleton._scream_if_not_simple()
            assert self._skeleton.is_connected(), "skeleton is not connected"
            assert self._skeleton.order() > 1, \
                "skeleton should have at least two vertices"
            edges = []
            blades = {frozenset(e): [] for e
                      in self._skeleton.edges(labels = False)}
            for f, l in faces.items():
                for i in range(len(l)):
                    e = frozenset([l[i-1], l[i]])
                    assert e in blades, "edge %s not in skeleton" % tuple(e)
                    if len(blades[e]) == 1:
                        g, = blades[e]
                        edges.append(((l[i-1], l[i], f), (l[i-1], l[i], g),
                                      FACE))
                        edges.append(((l[i], l[i-1], f), (l[i], l[i-1], g),
                                      FACE))
                    else:
                        assert len(blades[e]) == 0, \
                            "edge %s lies on three faces" % tuple(e)
                    blades[e].append(f)
                    edges.append(((l[i-1], l[i], f), (l[i], l[i-1], f),
                                  VERTEX))
                    edges.append(((l[i-1], l[i-2], f), (l[i-1], l[i], f),
                                  EDGE))
            assert all(len(x) == 2 for x in blades.values()), \
                "not all edges lie on two faces"
            kargs["data"] = edges
            vertex_fun = lambda (u, v, f): u
            edge_fun = lambda (u, v, f): frozenset([u, v])
            face_fun = lambda (u, v, f): f
            self._muscles = Graph(blades.values(), **NONSIMPLE)
        kargs["loops"] = False
        kargs["multiedges"] = False
        kargs["immutable"] = False
        edge_labels = kargs.pop("edge_labels", None)
        G = Graph(*largs, **kargs)
        assert G.is_regular(3), "graph is not 3-regular"
        if edge_labels is not None:
            for (u, v), l in edge_labels.items():
                G.set_edge_label(u, v, l)
        assert all(l in LABELS for l in G.edge_labels()), \
            "improper labels used"
        assert all(frozenset(G.edge_label(u, v) for v in G[u]) == LABELS
                   for u in G), "edges not labelled properly"
        edge_graph = Graph([(u, v) for u, v, l in G.edges() if l != EDGE])
        assert edge_graph.connected_components_number() * 4 == G.order(), \
            "involutions not commuting properly"
        Graph.__init__(self, G, immutable = True)
        if vertex_fun is None:
            self._vertex_fun = frozenset
        elif self._vertex_fun is None:
            self._vertex_fun = lambda v: vertex_fun(next(iter(v)))
        if edge_fun is None:
            self._edge_fun = frozenset
        elif self._edge_fun is None:
            self._edge_fun = lambda e: edge_fun(next(iter(e)))
        if face_fun is None:
            self._face_fun = frozenset
        elif self._face_fun is None:
            self._face_fun = lambda f: face_fun(next(iter(f)))
        self._vertices = {self._vertex_fun(x): frozenset(x) for x
                          in Graph([(u, v) for u, v, l in G.edges()
                                if l != VERTEX]).connected_components()}
        self._edges = {self._edge_fun(x): frozenset(x) for x
                       in edge_graph.connected_components()}
        self._faces = {self._face_fun(x): frozenset(x) for x
                       in Graph([(u, v) for u, v, l in G.edges()
                                 if l != FACE]).connected_components()}
        if self._skeleton is None:
            loops = {}
            dualloops = {}
            for e, t in self._edges.items():
                try:
                    loops[e] = next(v for v, s in self._vertices.items()
                                    if t <= s)
                except StopIteration:
                    pass
                try:
                    dualloops[e] = next(f for f, s in self._faces.items()
                                        if t <= s)
                except StopIteration:
                    pass
            self._skeleton = Graph([[loops[e]]*2 if e in loops
                                    else [v for v, s in self._vertices.items()
                                          if len(s & t) == 2]
                                    for e, t in self._edges.items()],
                                   **NONSIMPLE)
            self._muscles = Graph([[dualloops[e]]*2 if e in dualloops
                                   else [f for f, s in self._faces.items()
                                         if len(s & t) == 2]
                                   for e, t in self._edges.items()],
                                  **NONSIMPLE)
            if self._pos is not None:
                self._skeleton._pos = {}
                self._muscles._pos = {}
                for v in self._skeleton:
                    self._skeleton._pos[v] = meanpos(self, self._vertices[v])
                for f in self._muscles:
                    self._muscles._pos[f] = meanpos(self, self._faces[f])
        elif self._dual is not None:
            if self._dual._pos is not None:
                self._pos = dict(self._dual._pos)
        elif self._skeleton._pos is not None:
            if self._pos is None:
                self._pos = {}
                for u, v, f in self:
                    w = next(y for x, y, z in self[u, v, f]
                             if (x, z) == (u, f))
                    pu, pv, pw = (self._skeleton._pos[x] for x in (u, v, w))
                    self._pos[u, v, f] = tuple((1-C8-S8) * a +
                                               C8 * b + S8 * c
                                               for a, b, c in zip(pu, pv, pw))
            if self._muscles._pos is None:
                self._muscles._pos = {}
                for k, f in faces.items():
                    self._muscles._pos[k] = meanpos(self._skeleton, f)

    def characteristic(self):
        return self._skeleton.order() - self._skeleton.size() + \
            self._muscles.order()

    def copy(self, weighted = None, implementation = 'c_graph',
             data_structure = None, sparse = None, immutable = None):
        if immutable is False or (data_structure is not None
                                  and data_structure is not 'static_sparse'):
            return Graph(self).copy(weighted = weighted,
                                    implementation = implementation,
                                    data_structure = data_structure,
                                    sparse = sparse,
                                    immutable = immutable)
        else:
            return Graph.copy(self, weighted = weighted,
                                    implementation = implementation,
                                    data_structure = data_structure,
                                    sparse = sparse,
                                    immutable = immutable)

    def dual(self):
        if self._dual is None:
            self._dual = Tiling([(u, v, DUAL[l]) for u, v, l in self.edges()],
                                dual = self)
        return self._dual

    def graphplot(self, *largs, **kargs):
        if kargs.get("edge_colors") is None:
            kargs["edge_colors"] = {c: [(u, v) for u, v, l in self.edges()
                                        if l == t]
                                    for t, c in EDGE_COLORS.items()}
        return Graph.graphplot(self, *largs, **kargs)

    def muscles(self):
        return self._muscles

    def relabel(self, perm = None, inplace = True, return_map = False,
                check_input = True, complete_partial_function = True,
                immutable = True):
        if inplace:
            raise ValueError("To relabel an immutable graph use inplace=False")
        G = Graph(self, immutable = False)
        perm = G.relabel(perm, return_map = True, check_input = check_input,
                         complete_partial_function = complete_partial_function)
        if immutable is not False:
            G = self.__class__(self, vertex_labels = perm)
        if return_map:
            return G, perm
        else:
            return G

    def skeleton(self):
        return self._skeleton

    def _subgraph_by_adding(self, vertices = None, edges = None,
                            edge_property = None, immutable = None, *largs,
                            **kargs):
        if immutable is None:
            immutable = True
        return Graph(self)._subgraph_by_adding(vertices = vertices,
                                               edges = edges,
                                               edge_property = edge_property,
                                               immutable = immutable,
                                               *largs, **kargs)
