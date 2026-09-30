"""Exact rational model for the R1 finite response-network theorems.

One nonzero gain per oriented edge; reverse traversal uses the reciprocal.
No inner product, Hilbert completion, fitted constant, or floating point.
The supported graph is finite, connected, undirected, simple, and labelled.
"""

from collections import deque
from dataclasses import dataclass
from fractions import Fraction


def rational(value):
    """Accept exact inputs only; silently approximating floats is undesirable."""
    if not isinstance(value, (int, Fraction)):
        raise TypeError("Use integers or fractions.Fraction for exact inputs")
    return Fraction(value)


@dataclass(frozen=True)
class Edge:
    name: str
    tail: str
    head: str
    gain: Fraction

    def __post_init__(self):
        object.__setattr__(self, "gain", rational(self.gain))
        if self.gain == 0:
            raise ValueError("Reciprocal transport requires a nonzero gain")


class ResponseNetwork:
    def __init__(self, vertices, edges):
        self.vertices = tuple(vertices)
        self.edges = tuple(edges)
        if not self.vertices or len(set(self.vertices)) != len(self.vertices):
            raise ValueError("Vertices must be nonempty and distinct")
        if any(not isinstance(v, str) or not v for v in self.vertices):
            raise ValueError("Use nonempty string vertex labels")
        self._adj = {v: [] for v in self.vertices}
        self._edge = {}
        pairs = set()
        for edge in self.edges:
            if edge.name in self._edge:
                raise ValueError("Edge names must be distinct")
            if edge.tail not in self._adj or edge.head not in self._adj:
                raise ValueError("An edge endpoint is not a vertex")
            pair = frozenset((edge.tail, edge.head))
            if len(pair) != 2 or pair in pairs:
                raise ValueError("Self-loops and parallel edges are unsupported")
            pairs.add(pair)
            self._edge[edge.name] = edge
            self._adj[edge.tail].append((edge.head, edge.name, 1))
            self._adj[edge.head].append((edge.tail, edge.name, -1))
        self._tree(self.vertices[0], None)

    @property
    def cycle_rank(self):
        return len(self.edges) - len(self.vertices) + 1

    def _tree(self, root, tree_edges):
        if root not in self._adj:
            raise ValueError("Root is not a vertex")
        selected = None
        if tree_edges is not None:
            names = tuple(tree_edges)
            selected = set(names)
            if len(names) != len(selected) or len(names) != len(self.vertices) - 1:
                raise ValueError("A spanning tree needs n-1 distinct edges")
            if not selected.issubset(self._edge):
                raise ValueError("Unknown tree edge")
        parent = {root: None}
        order = []
        queue = deque((root,))
        while queue:
            here = queue.popleft()
            for there, name, sign in self._adj[here]:
                if selected is not None and name not in selected:
                    continue
                if there not in parent:
                    parent[there] = (here, name, sign)
                    order.append(there)
                    queue.append(there)
        if len(parent) != len(self.vertices):
            raise ValueError("Graph or chosen tree is disconnected")
        used = tuple(parent[v][1] for v in order)
        return parent, order, used

    def gain(self, tail, head):
        if tail not in self._adj or head not in self._adj:
            raise ValueError("Unknown vertex")
        for there, name, sign in self._adj[tail]:
            if there == head:
                gain = self._edge[name].gain
                return gain if sign == 1 else 1 / gain
        raise ValueError("No edge connects the requested vertices")

    def transport(self, walk):
        walk = tuple(walk)
        if not walk or any(v not in self._adj for v in walk):
            raise ValueError("A walk needs at least one known vertex")
        product = Fraction(1)
        for tail, head in zip(walk, walk[1:]):
            product *= self.gain(tail, head)
        return product

    def holonomy(self, walk):
        walk = tuple(walk)
        if not walk or walk[0] != walk[-1]:
            raise ValueError("A return walk must be closed")
        return self.transport(walk)

    def rescale(self, gauges):
        if set(gauges) != set(self.vertices):
            raise ValueError("Supply one nonzero gauge per vertex")
        gauges = {v: rational(gauges[v]) for v in self.vertices}
        if any(value == 0 for value in gauges.values()):
            raise ValueError("A reference rescaling must be invertible")
        return ResponseNetwork(self.vertices, (
            Edge(e.name, e.tail, e.head,
                 gauges[e.head] * e.gain / gauges[e.tail])
            for e in self.edges
        ))

    def normal_form(self, root=None, tree_edges=None):
        """Return (tree-normalized network, gauge factors, chord invariants)."""
        root = self.vertices[0] if root is None else root
        parent, order, used = self._tree(root, tree_edges)
        path_gain = {root: Fraction(1)}
        for vertex in order:
            previous, name, sign = parent[vertex]
            gain = self._edge[name].gain
            path_gain[vertex] = path_gain[previous] * gain ** sign
        gauges = {v: 1 / path_gain[v] for v in self.vertices}
        normalized = self.rescale(gauges)
        chords = {e.name: e.gain for e in normalized.edges if e.name not in used}
        return normalized, gauges, chords

    def path_operator(self, walk):
        """Matrix for the typed path action on the labelled direct sum Q^V."""
        walk = tuple(walk)
        product = self.transport(walk)
        index = {v: i for i, v in enumerate(self.vertices)}
        result = [[Fraction(0) for _ in self.vertices] for _ in self.vertices]
        result[index[walk[-1]]][index[walk[0]]] = product
        return result

    def rate_potential(self, rates, root=None, tree_edges=None):
        """Find eta with rate(i->j)=eta_j-eta_i, or return None.

        Rates are hypothetical instantaneous logarithmic gain derivatives.
        This verifies the algebraic condition, not an empirical time series.
        """
        if set(rates) != set(self._edge):
            raise ValueError("Supply one logarithmic rate per oriented edge")
        rates = {name: rational(rate) for name, rate in rates.items()}
        root = self.vertices[0] if root is None else root
        parent, order, _ = self._tree(root, tree_edges)
        potential = {root: Fraction(0)}
        for vertex in order:
            previous, name, sign = parent[vertex]
            potential[vertex] = potential[previous] + sign * rates[name]
        if any(rates[e.name] != potential[e.head] - potential[e.tail]
               for e in self.edges):
            return None
        return potential


def triangle():
    """Illustrative inputs, not constants derived from the framework."""
    return ResponseNetwork(('A', 'B', 'C'), (
        Edge('ab', 'A', 'B', 2),
        Edge('bc', 'B', 'C', 3),
        Edge('ca', 'C', 'A', Fraction(1, 7)),
    ))
