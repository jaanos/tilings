from sage.graphs.graph import Graph
from sage.misc.functional import numerical_approx as N
from sage.sets.set import Set
from .constants import VERTEX, EDGE, FACE, BLADE, CORNER
from .constants import LABELS, DUAL, EDGE_COLORS
from .constants import TRUNCATION_MAP, NONSIMPLE, T
from .functions import makeEdge, meanpos, flagPosition
from .functions import truncationVertexFunction

class Tiling(Graph):
    def __init__(self, *largs, **kargs):
        self._skeleton = None
        self._muscles = None
        self._vertex_fun = None
        self._edge_fun = None
        self._face_fun = None
        self._dual = kargs.pop("dual", None)
        vertex_rep = kargs.pop("vertex_rep", True)
        edge_rep = kargs.pop("edge_rep", True)
        face_rep = kargs.pop("face_rep", True)
        vertex_fun = kargs.pop("vertex_fun",
                               None if vertex_rep else Set)
        edge_fun = kargs.pop("edge_fun", None if edge_rep else Set)
        face_fun = kargs.pop("face_fun", None if face_rep else Set)
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
            blades = {Set(e): [] for e
                      in self._skeleton.edges(labels = False)}
            for f, l in faces.items():
                for i in range(len(l)):
                    e = Set([l[i-1], l[i]])
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
            edge_fun = lambda (u, v, f): Set([u, v])
            face_fun = lambda (u, v, f): f
            self._muscles = Graph(blades.values(), **NONSIMPLE)
        kargs["loops"] = False
        kargs["multiedges"] = True
        kargs["immutable"] = False
        edge_labels = kargs.pop("edge_labels", None)
        G = Graph(*largs, **kargs)
        assert G.is_regular(3), "graph is not 3-regular"
        if edge_labels is not None:
            for (u, v), l in edge_labels.items():
                G.set_edge_label(u, v, l)
        assert all(l in LABELS for l in G.edge_labels()), \
            "improper labels used"
        assert all(Set(sum((G.edge_label(u, v) for v in G[u]),
                            [])) == LABELS for u in G), \
            "edges not labelled properly"
        edge_graph = Graph([(u, v) for u, v, l in G.edges() if l != EDGE])
        assert edge_graph.connected_components_number() * 4 == G.order(), \
            "involutions not commuting properly"
        Graph.__init__(self, G, immutable = True, multiedges = True)
        if self._vertex_fun is None:
            if vertex_rep:
                self._vertex_fun = lambda v: vertex_fun(next(iter(v)))
            else:
                self._vertex_fun = vertex_fun
        if self._edge_fun is None:
            if edge_rep:
                self._edge_fun = lambda e: edge_fun(next(iter(e)))
            else:
                self._edge_fun = edge_fun
        if self._face_fun is None:
            if face_rep:
                self._face_fun = lambda f: face_fun(next(iter(f)))
            else:
                self._face_fun = face_fun
        self._vertices = {self._vertex_fun(x): Set(x) for x
                          in Graph([(u, v) for u, v, l in G.edges()
                                    if l != VERTEX],
                                   multiedges = True).connected_components()}
        self._edges = {self._edge_fun(x): Set(x) for x
                       in edge_graph.connected_components()}
        self._faces = {self._face_fun(x): Set(x) for x
                       in Graph([(u, v) for u, v, l in G.edges()
                                 if l != FACE],
                                 multiedges = True).connected_components()}
        if self._skeleton is None:
            self._skeleton = Graph([makeEdge([v for v, s
                                              in self._vertices.items()
                                              if len(s & t) > 0])
                                    for e, t in self._edges.items()],
                                   **NONSIMPLE)
            self._muscles = Graph([makeEdge([f for f, s in self._faces.items()
                                             if len(s & t) > 0])
                                   for e, t in self._edges.items()],
                                  **NONSIMPLE)
            if self._pos is not None:
                self._reposition()
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
                    self._pos[u, v, f] = flagPosition(pu, pv, pw)
            if self._muscles._pos is None:
                self._muscles._pos = {}
                for k, f in faces.items():
                    self._muscles._pos[k] = meanpos(self._skeleton, f)

    def automorphism_group(self, partition = None, verbosity = 0,
                           edge_labels = True, order = False,
                           return_group = True, orbits = False,
                           algorithm = None):
        return Graph.automorphism_group(self, partition = partition,
                                        edge_labels = edge_labels,
                                        order = order,
                                        return_group = return_group,
                                        orbits = orbits,
                                        algorithm = algorithm)

    def canonical_label(self, partition = None, certificate = False,
                        edge_labels = True, algorithm = None,
                        return_graph = True):
        return Graph.canonical_label(self, partition = partition,
                                     certificate = certificate,
                                     edge_labels = edge_labels,
                                     algorithm = algorithm,
                                     return_graph = return_graph)

    def chamfer(self):
        edges = []
        for s in self:
            edges.append(((s, FACE), (s, EDGE), FACE))
            edges.append(((s, CORNER), (s, EDGE), EDGE))
            edges.append(((s, CORNER), (s, VERTEX), VERTEX))
        for s, t, l in self.edges(labels = True):
            if l == VERTEX:
                edges.append(((s, EDGE), (t, EDGE), VERTEX))
                edges.append(((s, FACE), (t, FACE), VERTEX))
            elif l == EDGE:
                edges.append(((s, VERTEX), (t, VERTEX), FACE))
                edges.append(((s, CORNER), (t, CORNER), FACE))
                edges.append(((s, FACE), (t, FACE), EDGE))
            elif l == FACE:
                edges.append(((s, VERTEX), (t, VERTEX), EDGE))
        return Tiling(edges, pos = self._chamferPosition(),
                      vertex_rep = False, edge_rep = False, face_rep = False,
                      vertex_fun = self._chamferVertexFunction,
                      edge_fun = self._chamferEdgeFunction,
                      face_fun = self._chamferFaceFunction)

    def _chamferVertexFunction(self, s):
        p, l = next(iter(s))
        r = truncationVertexFunction(s)
        if l == VERTEX:
            return (self._vertex_fun(r), VERTEX)
        else:
            return (r, CORNER)

    def _chamferEdgeFunction(self, s):
        p, l = next(iter(s))
        r = truncationVertexFunction(s)
        if l in [VERTEX, CORNER]:
            return (r, CORNER)
        else:
            return (r, BLADE)

    def _chamferFaceFunction(self, s):
        p, l = next(iter(s))
        r = truncationVertexFunction(s)
        if l == FACE:
            return (self._face_fun(r), FACE)
        else:
            return (r, EDGE)

    def _chamferPosition(self):
        if self._pos is None:
            return None
        vertex = {}
        neigh = {}
        vpos = {}
        pos = {(s, VERTEX): p for s, p in self._pos.items()}
        for v, c in self._vertices.items():
            for t in c:
                vertex[t] = v
        for s, t, l in self.edges(labels = True):
            for r in [s, t]:
                if r not in neigh:
                    neigh[r] = {}
            neigh[s][l] = t
            neigh[t][l] = s
            if l == EDGE:
                x = [p+q-r for p, q, r in zip(self._pos[s], self._pos[t],
                                              self._skeleton._pos[vertex[s]])]
                vpos[s] = x
                vpos[t] = x
        for s in self:
            n = neigh[s]
            p = vpos[s]
            q = vpos[neigh[s][VERTEX]]
            r = self._skeleton._pos[vertex[s]]
            pos[s, CORNER] = flagPosition(p, r, q)
            pos[s, EDGE] = flagPosition(p, q, r)
            pos[s, FACE] = flagPosition(p, q,
                                        vpos[neigh[neigh[s][EDGE]][VERTEX]])
        return pos

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
            kargs["edge_colors"] = {c: [(u, v, l) for u, v, l in self.edges()
                                        if l == t]
                                    for t, c in EDGE_COLORS.items()}
        return Graph.graphplot(self, *largs, **kargs)

    def is_isomorphic(self, other, certificate = False, edge_labels = True):
        return Graph.is_isomorphic(self, other, certificate = certificate,
                                   edge_labels = edge_labels)

    def is_vertex_transitive(self, partition = None, edge_labels = True,
                             order = False, return_group = True,
                             orbits = False):
        return Graph.is_vertex_transitive(self, partition = partition,
                                          edge_labels = edge_labels,
                                          order = order,
                                          return_group = return_group,
                                          orbits = orbits)

    def muscles(self):
        return self._muscles

    def relabel(self, perm = None, inplace = True, return_map = False,
                check_input = True, complete_partial_function = True,
                immutable = True):
        if inplace:
            raise ValueError("To relabel an immutable graph use inplace=False")
        G = Graph(self, immutable = False)
        return G.relabel(perm, inplace = inplace, return_map = return_map,
                         check_input = check_input, immutable = immutable,
                         complete_partial_function = complete_partial_function)

    def _reposition(self):
        self._skeleton._pos = {}
        self._muscles._pos = {}
        for v in self._skeleton:
            self._skeleton._pos[v] = meanpos(self, self._vertices[v])
        for f in self._muscles:
            self._muscles._pos[f] = meanpos(self, self._faces[f])

    def set_pos(self, pos, dim = 2):
        Graph.set_pos(self, pos, dim = 2)
        self._reposition()

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

    def truncation(self):
        edges = []
        pos = None
        for s in self:
            edges.append(((s, EDGE), (s, FACE), EDGE))
            edges.append(((s, VERTEX), (s, FACE), FACE))
        for s, t, l in self.edges(labels = True):
            if l == VERTEX:
                edges.append(((s, EDGE), (t, EDGE), VERTEX))
            elif l == EDGE:
                edges.append(((s, FACE), (t, FACE), VERTEX))
                edges.append(((s, VERTEX), (t, VERTEX), VERTEX))
            elif l == FACE:
                edges.append(((s, EDGE), (t, EDGE), FACE))
                edges.append(((s, VERTEX), (t, VERTEX), EDGE))
        if self._pos is not None:
            pos = {}
            for s, t, l in self.edges(labels = True):
                if l != FACE:
                    for a, b in [(s, t), (t, s)]:
                        pos[a, TRUNCATION_MAP[l]] = [N((1-T)*p + T*q)
                                                     for p, q in
                                                     zip(self._pos[a],
                                                         self._pos[b])]
            for s in self:
                pos[s, FACE] = [p+q-r for p, q, r in
                                zip(pos[s, VERTEX], pos[s, EDGE],
                                    self._pos[s])]
        return Tiling(edges, pos = pos,
                      vertex_rep = False, edge_rep = False, face_rep = False,
                      vertex_fun = truncationVertexFunction,
                      edge_fun = self._truncationEdgeFunction,
                      face_fun = self._truncationFaceFunction)

    def _truncationEdgeFunction(self, s):
        p, l = next(iter(s))
        r = truncationVertexFunction(s)
        if l == EDGE:
            return (self._edge_fun(r), EDGE)
        else:
            return (r, CORNER)

    def _truncationFaceFunction(self, s):
        p, l = next(iter(s))
        r = truncationVertexFunction(s)
        if l == VERTEX:
            return (self._vertex_fun(r), VERTEX)
        else:
            return (self._face_fun(r), FACE)
