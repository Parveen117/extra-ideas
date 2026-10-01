"""R10: sourced native Bianchi and an explicitly admitted Riemann adapter.

The metric adapter is a boundary contract, not a primitive of the native
algebra. All arithmetic is rational; source geometry remains unchanged.
Connection jets, cut grading and finite transports reuse R7-R9.
"""

from dataclasses import dataclass
from fractions import Fraction as Q
from aghora_return import inverse, transpose
from emk_curvature_balance import grade
from emk_curvature_observation import (
    ConnectionJet, K, R, RK, add, commutator, identity, is_zero,
    matrix, mul, same_carrier, scalar, scale, sub,
)
from emk_tensor_calculus import Move, nullspace


def rectangular(rows, n, m):
    value = tuple(tuple(scalar(x) for x in row) for row in rows)
    if len(value) != n or any(len(row) != m for row in value):
        raise ValueError('Rectangular map dimensions do not match the declared carriers')
    return value


def rzero(n, m):
    return tuple(tuple(Q(0) for _ in range(m)) for _ in range(n))


def rproduct(a, b):
    if not a or not b or len(a[0]) != len(b):
        raise ValueError('Rectangular composition has no lawful matching middle carrier')
    return tuple(tuple(sum((a[i][k]*b[k][j] for k in range(len(b))), Q(0))
                       for j in range(len(b[0]))) for i in range(len(a)))


def rsubtract(a, b):
    if len(a) != len(b) or any(len(ar) != len(br) for ar, br in zip(a, b)):
        raise ValueError('Rectangular subtraction requires equal types')
    return tuple(tuple(x-y for x, y in zip(ar, br)) for ar, br in zip(a, b))


def sourced_bianchi(operators, cut):
    """Constant-coefficient cyclic law J_even + J_hidden = 0.

    The smooth version also has derivatives. This exact finite function
    deliberately declares constant coefficients and a constant cut.
    """
    if len(operators) != 3:
        raise ValueError('Declare exactly three comparison directions')
    a = same_carrier(*operators)
    even = tuple(grade(x, cut)[0] for x in a)
    odd = tuple(grade(x, cut)[1] for x in a)
    raw = scale(a[0], 0)
    visible, hidden = raw, raw
    curvatures = []
    for i, j, k in ((0, 1, 2), (1, 2, 0), (2, 0, 1)):
        f = commutator(a[j], a[k])
        fe, fo = grade(f, cut)
        curvatures.append(f)
        raw = add(raw, commutator(a[i], f))
        visible = add(visible, commutator(even[i], fe))
        hidden = add(hidden, commutator(odd[i], fo))
    return {'cyclic_curvatures': tuple(curvatures), 'full_bianchi': raw,
            'visible_bianchi': visible, 'hidden_source': hidden,
            'sourced_residual': add(visible, hidden)}


def sourced_bianchi_jet(connection, second_partials, cut):
    """Three-direction smooth Bianchi at a point, for a constant cut.

    second_partials[i][j][k]=d_k d_j A_i. Commuting coordinate partials
    are a checked contract; no derivative of a moving cut is omitted.
    """
    if len(connection.directions) != 3:
        raise ValueError('The cyclic jet certificate requires three comparison directions')
    if (len(second_partials) != 3 or any(len(row) != 3 for row in second_partials)
            or any(len(row) != 3 for plane in second_partials for row in plane)):
        raise ValueError('Declare all second connection partials')
    second = tuple(tuple(tuple(matrix(x) for x in row) for row in plane) for plane in second_partials)
    same_carrier(*connection.operators, *(x for plane in second for row in plane for x in row))
    if any(second[i][j][k] != second[i][k][j] for i in range(3) for j in range(3) for k in range(3)):
        raise ValueError('Coordinate second partials must commute')
    even, odd = zip(*(grade(a, cut) for a in connection.operators))
    raw = scale(connection.operators[0], 0)
    visible, hidden = raw, raw
    derivatives = []
    for i, j, k in ((0, 1, 2), (1, 2, 0), (2, 0, 1)):
        f = connection.curvature(j, k)
        df = add(sub(second[k][j][i], second[j][k][i]),
                 add(commutator(connection.partials[j][i], connection.operators[k]),
                     commutator(connection.operators[j], connection.partials[k][i])))
        fe, fo = grade(f, cut)
        derivatives.append(df)
        raw = add(raw, add(df, commutator(connection.operators[i], f)))
        visible = add(visible, add(grade(df, cut)[0], commutator(even[i], fe)))
        hidden = add(hidden, commutator(odd[i], fo))
    return {'cyclic_curvature_derivatives': tuple(derivatives), 'full_bianchi': raw,
            'visible_bianchi': visible, 'hidden_source': hidden,
            'sourced_residual': add(visible, hidden)}


@dataclass(frozen=True)
class MetricJet:
    """Admitted nondegenerate symmetric metric and first/second coordinate jets.

    first[i]=d_i g; second[i][j]=d_j d_i g. Both metric and derivative
    symmetries are validated; this is a local, not a global PDE certificate.
    """
    directions: tuple
    value: tuple
    first: tuple
    second: tuple

    def __post_init__(self):
        names, g = tuple(self.directions), matrix(self.value)
        n = len(g)
        if (len(names) != n or len(set(names)) != n
                or any(not isinstance(x, str) or not x for x in names)):
            raise ValueError('The Riemann adapter needs one tangent mode per coordinate direction')
        if transpose(g) != g:
            raise ValueError('The admitted metric must be symmetric')
        inverse(g)
        if len(self.first) != n or len(self.second) != n or any(len(row) != n for row in self.second):
            raise ValueError('Metric derivative jet dimensions are inconsistent')
        first = tuple(matrix(x) for x in self.first)
        second = tuple(tuple(matrix(x) for x in row) for row in self.second)
        same_carrier(g, *first, *(x for row in second for x in row))
        if (any(transpose(x) != x for x in first)
                or any(transpose(second[i][j]) != second[i][j] or second[i][j] != second[j][i]
                       for i in range(n) for j in range(n))):
            raise ValueError('A metric two-jet must have symmetric matrices and commuting partials')
        object.__setattr__(self, 'directions', names)
        object.__setattr__(self, 'value', g)
        object.__setattr__(self, 'first', first)
        object.__setattr__(self, 'second', second)

    def levi_civita(self):
        n, g = len(self.value), self.value
        gi = inverse(g)
        derivative_inverse = tuple(scale(mul(mul(gi, dg), gi), -1) for dg in self.first)

        def coefficient(i, a, b):
            return sum((gi[a][c]*(self.first[i][c][b]+self.first[b][c][i]-self.first[c][i][b])
                        for c in range(n)), Q(0))/2

        def derivative(i, j, a, b):
            return sum((
                derivative_inverse[j][a][c]*(self.first[i][c][b]+self.first[b][c][i]-self.first[c][i][b])
                + gi[a][c]*(self.second[i][j][c][b]+self.second[b][j][c][i]-self.second[c][j][i][b])
                for c in range(n)), Q(0))/2

        operators = tuple(matrix([[coefficient(i, a, b) for b in range(n)] for a in range(n)]) for i in range(n))
        partials = tuple(tuple(matrix([[derivative(i, j, a, b) for b in range(n)] for a in range(n)])
                               for j in range(n)) for i in range(n))
        return ConnectionJet(self.directions, operators, partials)


def warp_metric(kappa, v):
    """Reuse EMK-G1's declared W=1+kappa*v^2 with W>0; no square root."""
    kappa, v = scalar(kappa), scalar(v)
    w, wp, wpp = 1+kappa*v*v, 2*kappa*v, 2*kappa
    if w <= 0:
        raise ValueError('The selected seam metric domain requires W>0')
    z = scale(identity(2), 0)
    return MetricJet(('u', 'v'), matrix([[w, 0], [0, 1]]),
                     (z, matrix([[wp, 0], [0, 0]])),
                     ((z, z), (z, matrix([[wpp, 0], [0, 0]]))))


def flat_metric(n=2, directions=None):
    names = tuple(directions) if directions is not None else tuple(f'x{i}' for i in range(n))
    z = scale(identity(n), 0)
    return MetricJet(names, identity(n), (z,)*n, tuple((z,)*n for _ in range(n)))


def tangent_compatibility(connection, metric):
    n = len(metric.value)
    if connection.directions != metric.directions or len(connection.operators[0]) != n:
        raise ValueError('The tangent connection and metric must have the same coordinate contract')
    nonmetricity = tuple(sub(metric.first[i], add(mul(transpose(a), metric.value), mul(metric.value, a)))
                        for i, a in enumerate(connection.operators))
    torsion = tuple(tuple(tuple(connection.operators[i][a][j]-connection.operators[j][a][i]
                                 for a in range(n)) for j in range(n)) for i in range(n))
    lc = metric.levi_civita()
    values = tuple(sub(a, g) for a, g in zip(connection.operators, lc.operators))
    partials = tuple(tuple(sub(a, g) for a, g in zip(ar, gr)) for ar, gr in zip(connection.partials, lc.partials))
    return {'nonmetricity': nonmetricity, 'torsion': torsion,
            'distortion_values': values, 'distortion_partials': partials,
            'Levi_Civita_connection_jet': all(is_zero(x) for x in values) and all(is_zero(x) for row in partials for x in row)}


def projection(n, total):
    if not 0 < n <= total:
        raise ValueError('Visible tangent dimension must fit inside the native carrier')
    return rectangular([[int(i == j) for j in range(total)] for i in range(n)], n, total)


def descent_report(native, visible, observer):
    """Constant, surjective C; test both value AND differentiated intertwining."""
    if native.directions != visible.directions:
        raise ValueError('Connection comparison directions do not match')
    total, n = len(native.operators[0]), len(visible.operators[0])
    c = rectangular(observer, n, total)
    if len(nullspace(c, total)) != total-n:
        raise ValueError('A tangent observer must be surjective')
    values = tuple(rsubtract(rproduct(c, a), rproduct(g, c)) for a, g in zip(native.operators, visible.operators))
    partials = tuple(tuple(rsubtract(rproduct(c, a), rproduct(g, c)) for a, g in zip(ar, gr))
                     for ar, gr in zip(native.partials, visible.partials))
    curvatures = tuple(rsubtract(rproduct(c, native.curvature(i, j)), rproduct(visible.curvature(i, j), c))
                       for i in range(len(native.directions)) for j in range(i+1, len(native.directions)))
    return {'value_intertwining': values, 'derivative_intertwining': partials,
            'curvature_intertwining': curvatures,
            'connection_jet_descends': all(is_zero(x) for x in values) and all(is_zero(x) for row in partials for x in row),
            'curvature_descends': all(is_zero(x) for x in curvatures)}


def riemann_reduction(native, metric, observer):
    lc = metric.levi_civita()
    report = descent_report(native, lc, observer)
    report['Riemann_sector_certified'] = report['connection_jet_descends'] and report['curvature_descends']
    report['tangent_curvature'] = lc.components()
    return report


def block(a, upper, lower, h):
    a, h = matrix(a), matrix(h)
    n, m = len(a), len(h)
    upper, lower = rectangular(upper, n, m), rectangular(lower, m, n)
    return matrix([a[i]+upper[i] for i in range(n)]+[lower[i]+h[i] for i in range(m)])


def extend_connection(visible, hidden, upper=None, lower=None):
    """Admitted block extension; off-diagonal couplings are constant here."""
    if visible.directions != hidden.directions:
        raise ValueError('Both sectors must share base comparison directions')
    n, m, d = len(visible.operators[0]), len(hidden.operators[0]), len(visible.directions)
    upper = (rzero(n, m),)*d if upper is None else tuple(rectangular(x, n, m) for x in upper)
    lower = (rzero(m, n),)*d if lower is None else tuple(rectangular(x, m, n) for x in lower)
    if len(upper) != d or len(lower) != d:
        raise ValueError('Declare one upper and lower coupling per direction')
    operators = tuple(block(a, e, l, h) for a, e, l, h in zip(visible.operators, upper, lower, hidden.operators))
    partials = tuple(tuple(block(visible.partials[i][j], rzero(n, m), rzero(m, n), hidden.partials[i][j])
                           for j in range(d)) for i in range(d))
    return ConnectionJet(visible.directions, operators, partials)


def visible_decomposition(native, metric, i=0, j=1):
    """Visible full curvature = Riemann + distortion + hidden excursions."""
    if native.directions != metric.directions:
        raise ValueError('Native and metric direction declarations differ')
    n, total = len(metric.value), len(native.operators[0])
    if total < n:
        raise ValueError('Native carrier is smaller than the tangent target')
    corner = lambda a: matrix([row[:n] for row in a[:n]])
    visible = ConnectionJet(native.directions, tuple(map(corner, native.operators)),
                            tuple(tuple(map(corner, row)) for row in native.partials))
    lc = metric.levi_civita()
    s = tuple(sub(a, g) for a, g in zip(visible.operators, lc.operators))
    ds = tuple(tuple(sub(a, g) for a, g in zip(ar, gr)) for ar, gr in zip(visible.partials, lc.partials))
    distortion = add(sub(ds[j][i], ds[i][j]),
                     add(add(commutator(lc.operators[i], s[j]), commutator(s[i], lc.operators[j])),
                         commutator(s[i], s[j])))
    excursion = scale(metric.value, 0)
    if total > n:
        e_i = tuple(row[n:] for row in native.operators[i][:n])
        e_j = tuple(row[n:] for row in native.operators[j][:n])
        l_i = tuple(row[:n] for row in native.operators[i][n:])
        l_j = tuple(row[:n] for row in native.operators[j][n:])
        excursion = matrix(rsubtract(rproduct(e_i, l_j), rproduct(e_j, l_i)))
    riemann = lc.curvature(i, j)
    observed = corner(native.curvature(i, j))
    return {'visible_full_curvature': observed, 'Riemann': riemann,
            'distortion': distortion, 'excursion': excursion,
            'reconstruction_residual': sub(observed, add(riemann, add(distortion, excursion)))}


def lifted_order_loop():
    """Finite transport witness, separate from infinitesimal connection curvature."""
    zero = rzero(2, 2)
    ur, uk = block(identity(2), zero, zero, R), block(identity(2), zero, zero, K)
    path = Move('0', 'a', (('V', ur),)).then(Move('a', '2', (('V', uk),)))
    alternate = Move('0', 'b', (('V', uk),)).then(Move('b', '2', (('V', ur),)))
    return path.then(alternate.reverse())


def exact_witnesses(geometry):
    def elementary(i, j):
        return matrix([[int(a == i and bb == j) for bb in range(3)] for a in range(3)])
    cut = matrix([[1, 0, 0], [0, 1, 0], [0, 0, -1]])
    sourced = sourced_bianchi((elementary(0, 1), elementary(1, 2), elementary(2, 0)), cut)
    metric = warp_metric(Q(3), Q(1, 3))
    lc = metric.levi_civita()
    hidden = ConnectionJet.constant(metric.directions, (K, R))
    native = extend_connection(lc, hidden, lower=(K, R))
    c = projection(2, 4)
    flat = flat_metric()
    hidden_flat = ConnectionJet.constant(flat.directions, (K, R))
    flat_lift = extend_connection(flat.levi_civita(), hidden_flat)
    feedback = extend_connection(flat.levi_civita(), hidden_flat, upper=(identity(2), scale(K, 0)),
                                 lower=(scale(K, 0), K))
    z = scale(identity(2), 0)
    point_only = ConnectionJet(flat.directions, (z, z), ((z, K), (z, z)))
    loop = lifted_order_loop()
    return {'sourced_Bianchi': sourced,
            'Riemann_descent': riemann_reduction(native, metric, c),
            'curved_metric_visible_decomposition': visible_decomposition(native, metric),
            'source_Gaussian_curvature': geometry.gaussian_curvature_closed(Q(1, 3), Q(3)),
            'flat_Riemann_with_curved_hidden_sector': {
                'tangent_curvature': flat.levi_civita().curvature(0, 1),
                'full_native_curvature': flat_lift.curvature(0, 1),
                'reduction': riemann_reduction(flat_lift, flat, c)},
            'feedback_breaks_descent': {'decomposition': visible_decomposition(feedback, flat),
                                       'descent': riemann_reduction(feedback, flat, c)},
            'pointwise_metric_torsion_test_is_insufficient': {
                'compatibility': tangent_compatibility(point_only, flat),
                'actual_curvature': point_only.curvature(0, 1),
                'reduction': riemann_reduction(point_only, flat, identity(2))},
            'finite_loop': {'full_carrier': dict(loop.operators)['V'],
                            'visible_carrier': matrix([row[:2] for row in dict(loop.operators)['V'][:2]])}}
