from sage.graphs.graph import Graph

LABELS = frozenset(["v", "e", "f"])
EDGE_COLORS = {"v": "red", "e": "green", "f": "blue"}

class Tiling(Graph):
    def __init__(self, *largs, **kargs):
        kargs["loops"] = False
        kargs["multiedges"] = False
        kargs["immutable"] = False
        edge_labels = None
        if "edge_labels" in kargs:
            edge_labels = kargs["edge_labels"]
            del kargs["edge_labels"]
        G = Graph(*largs, **kargs)
        assert G.is_regular(3), "graph is not 3-regular"
        if edge_labels is not None:
            for (u, v), l in edge_labels.items():
                G.set_edge_label(u, v, l)
        assert all(l in LABELS for l in G.edge_labels()), \
            "improper labels used"
        assert all(frozenset(G.edge_label(u, v) for v in G[u]) == LABELS
                   for u in G), "edges not labelled properly"
        assert Graph([(u, v) for u, v, l in G.edges()
                      if l != "e"]).connected_components_number() * 4 == \
                G.order(), "involutions not commuting properly"
        Graph.__init__(self, G, immutable = True)

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

    def graphplot(self, *largs, **kargs):
        if kargs.get("edge_colors") is None:
            kargs["edge_colors"] = {c: [(u, v) for u, v, l in self.edges()
                                        if l == t]
                                    for t, c in EDGE_COLORS.items()}
        return Graph.graphplot(self, *largs, **kargs)

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
