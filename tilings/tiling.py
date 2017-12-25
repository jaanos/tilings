from sage.graphs.graph import Graph

VERTEX = "v"
EDGE = "e"
FACE = "f"
LABELS = frozenset([VERTEX, EDGE, FACE])
DUAL = {VERTEX: FACE, EDGE: EDGE, FACE: VERTEX}
EDGE_COLORS = {VERTEX: "red", EDGE: "green", FACE: "blue"}
NONSIMPLE = {"loops": True, "multiedges": True}

class Tiling(Graph):
    def __init__(self, *largs, **kargs):
        self._skeleton = None
        self._muscles = None
        self._dual = kargs.pop("dual", None)
        if self._dual is not None:
            self._skeleton = self._dual._muscles
            self._muscles = self._dual._skeleton
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
            self._vertices = self._skeleton.vertices()
            self._edges = self._skeleton.edges(labels = False)
            self._faces = faces.keys()
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
        if self._skeleton is None:
            self._vertices = [frozenset(x) for x
                              in Graph([(u, v) for u, v, l in G.edges()
                                    if l != VERTEX]).connected_components()]
            self._edges = [frozenset(x) for x
                           in edge_graph.connected_components()]
            self._faces = [frozenset(x) for x
                           in Graph([(u, v) for u, v, l in G.edges()
                                     if l != FACE]).connected_components()]
            loops = {}
            dualloops = {}
            for e in self._edges:
                try:
                    loops[e] = next(v for v in self._vertices
                                    if e.issubset(v))
                except StopIteration:
                    pass
                try:
                    dualloops[e] = next(v for v in self._faces
                                        if e.issubset(v))
                except StopIteration:
                    pass
            self._skeleton = Graph([[loops[e], loops[e]] if e in loops
                                    else [v for v in self._vertices
                                          if len(e & v) == 2]
                                    for e in self._edges], **NONSIMPLE)
            self._muscles = Graph([[dualloops[e], dualloops[e]]
                                   if e in dualloops
                                   else [f for f in self._faces
                                         if len(e & f) == 2]
                                   for e in self._edges], **NONSIMPLE)

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
