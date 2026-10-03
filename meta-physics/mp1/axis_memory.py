"""MP-1: exact axis-observer memory on the existing rational EMK carrier.

No floating point, fitted dynamics, external runtime, or new native algebra.
The continuous claims and the between-sample contracts are in THEOREM.md.
"""

from fractions import Fraction as Q
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / '04-operator-evolution'))
from emk_tensor_calculus import I, R, K, RK  # noqa: E402
from aghora_return import add, scale, mul, inverse  # noqa: E402


def rational(x):
    if type(x) is not int and not isinstance(x, Q):
        raise TypeError('Use exact integers or Fraction values, not floats/bools')
    return Q(x)


def point(p):
    p = tuple(rational(x) for x in p)
    if len(p) != 2:
        raise ValueError('A coefficient pair is required')
    return p


def dot(a, b):
    return a[0] * b[0] + a[1] * b[1]


def det(a, b):
    return a[0] * b[1] - a[1] * b[0]


def unit_pair(p):
    p = point(p)
    if dot(p, p) != 1:
        raise ValueError('A unit native turn is required')
    return p


def rotor(p):
    a, b = unit_pair(p)
    return add(scale(I, a), scale(R, b))


def rational_turn(t):
    """Rational unit turn, with (-1, 0) supplied separately at chart infinity."""
    t = rational(t)
    return ((1 - t*t) / (1 + t*t), 2*t / (1 + t*t))


def axis_coefficients(p):
    a, b = unit_pair(p)
    return (a*a - b*b, 2*a*b)


def axis_operator(p):
    g = rotor(p)
    return mul(mul(g, K), inverse(g))


def reference_intensity(p, v=(1, 0)):
    """Declared retained-reference probe ||v+g v||^2; no Born-law premise."""
    v = point(v)
    g = rotor(p)
    out = tuple(v[i] + sum(g[i][j]*v[j] for j in range(2)) for i in range(2))
    return dot(out, out)


def segment_clearance_squared(a, b):
    a, b = point(a), point(b)
    d = (b[0] - a[0], b[1] - a[1])
    dd = dot(d, d)
    if dd == 0:
        return dot(a, a)
    t = min(Q(1), max(Q(0), -dot(a, d) / dd))
    c = (a[0] + t*d[0], a[1] + t*d[1])
    return dot(c, c)


def closed_polygon(vertices):
    p = tuple(point(v) for v in vertices)
    if len(p) < 2 or p[0] != p[-1]:
        raise ValueError('An explicitly closed polygon is required')
    if any(segment_clearance_squared(a, b) == 0 for a, b in zip(p, p[1:])):
        raise ValueError('A vertex or segment meets zero; phase is undefined')
    return p


def winding(vertices):
    """Oriented positive-ray crossing count with half-open vertex convention."""
    p = closed_polygon(vertices)
    count = 0
    for a, b in zip(p, p[1:]):
        cross = det(a, b)
        if a[1] <= 0 < b[1] and cross > 0:
            count += 1
        elif b[1] <= 0 < a[1] and cross < 0:
            count -= 1
    return count


def polygon_report(vertices, *, vertex_error=Q(0), interpolation_error=Q(0)):
    """Certify a polygon and a DECLARED coordinatewise error tube.

    The true path must be closed. Its vertices are within vertex_error of
    the measured vertices; between samples it is within interpolation_error
    of its own straight chord, coordinatewise. These are caller hypotheses,
    not facts inferred from samples. Without them only the polygon is certified.
    """
    p = closed_polygon(vertices)
    eps, eta = rational(vertex_error), rational(interpolation_error)
    if eps < 0 or eta < 0:
        raise ValueError('Error budgets must be nonnegative')
    clearance = min(segment_clearance_squared(a, b) for a, b in zip(p, p[1:]))
    tube_squared = 2 * (eps + eta)**2
    nu = winding(p)
    return {
        'winding': nu,
        'relative_frame_sign': -1 if nu % 2 else 1,
        'clearance_squared': clearance,
        'tube_squared_bound': tube_squared,
        'tube_certified': tube_squared < clearance,
        'scope': 'closed rational polygon; continuous path only under declared tube contract',
    }


def principal_frame_lift(axis_turns, initial_frame_turn=Q(0)):
    """Conditional sample recovery, in full-turn units.

    Valid for the actual frame only if each UNWRAPPED frame increment has
    absolute value <1/4. Samples cannot establish that hypothesis. Exactly
    opposite axis samples are rejected, never rounded to a guessed sign.
    """
    axis = tuple(rational(x) for x in axis_turns)
    if not axis or any(not 0 <= x < 1 for x in axis):
        raise ValueError('Axis samples must be wrapped turns in [0,1)')
    start = rational(initial_frame_turn)
    if (2*start) % 1 != axis[0]:
        raise ValueError('Initial frame and axis do not match')
    out = [start]
    for a, b in zip(axis, axis[1:]):
        delta = (b - a) % 1
        if delta == Q(1, 2):
            raise ValueError('Half-turn axis tie: frame sheet is ambiguous')
        if delta > Q(1, 2):
            delta -= 1
        out.append(out[-1] + delta / 2)
    return tuple(out)


def required_history_labels(max_abs_winding):
    if type(max_abs_winding) is not int or max_abs_winding < 0:
        raise ValueError('A nonnegative integer winding budget is required')
    return 2*max_abs_winding + 1


def diamond(turns=1):
    """Exact example, also used to exhibit hidden endpoint aliases."""
    if type(turns) is not int:
        raise TypeError('An integer turn count is required')
    cycle = ((1, 0), (0, 1), (-1, 0), (0, -1), (1, 0))
    if turns < 0:
        cycle = tuple(reversed(cycle))
    return cycle[:-1] * abs(turns) + (cycle[0],) if turns else ((1, 0), (1, 0))
