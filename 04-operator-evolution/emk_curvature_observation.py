"""R8 research adapter: exact curvature, compression and observation.

Uses the existing R2 rational matrix implementation and R7 KIR carrier.
The canonical native algebra engine remains in RKF/operator_foundation.
ConnectionJet is an admitted coordinate adapter at ONE point, not a PDE
solver or a universal definition of the full EMK master tensor.
"""

from dataclasses import dataclass
from fractions import Fraction as Q

from aghora_return import (add, conjugate, identity, inverse, matrix, mul,
                          parity_projectors, scalar, scale, sub)
from emk_tensor_calculus import K, R, RK


def same_carrier(*operators):
    result = tuple(matrix(a) for a in operators)
    if not result or len({len(a) for a in result}) != 1:
        raise ValueError('Operators must act on one nonempty square carrier')
    return result


def commutator(a, b):
    a, b = same_carrier(a, b)
    return sub(mul(a, b), mul(b, a))


def is_zero(a):
    return all(x == 0 for row in a for x in row)


def trace(a):
    a = matrix(a)
    return sum((a[i][i] for i in range(len(a))), Q(0))


def compression_report(a, b, projector):
    """C=[PAP,PBP]=P[A,B]P-E; E=PAQBP-PBQAP, Q=I-P.

    Idempotence is sufficient. No positivity or orthogonality is assumed.
    Compression is an intervention between compositions, not a homomorphism.
    """
    a, b, p = same_carrier(a, b, projector)
    if mul(p, p) != p:
        raise ValueError('Compression requires an idempotent projector')
    q = sub(identity(len(p)), p)
    full = commutator(a, b)
    readout = mul(mul(p, full), p)
    reduced = commutator(mul(mul(p, a), p), mul(mul(p, b), p))
    excursions = sub(mul(mul(mul(mul(p, a), q), b), p),
                     mul(mul(mul(mul(p, b), q), a), p))
    return {'full': full, 'compressed_full': readout,
            'reduced': reduced, 'excursions': excursions,
            'identity_residual': sub(reduced, sub(readout, excursions)),
            'full_flat': is_zero(full), 'reduced_flat': is_zero(reduced)}


def pair(covector, operator, vector):
    operator = matrix(operator)
    left, right = tuple(map(scalar, covector)), tuple(map(scalar, vector))
    if len(left) != len(operator) or len(right) != len(operator):
        raise ValueError('Covector, carrier and vector dimensions must agree')
    return sum((left[i] * operator[i][j] * right[j]
                for i in range(len(operator)) for j in range(len(operator))), Q(0))


W_PLUS, W_MINUS = (Q(1), Q(1)), (Q(1), Q(-1))
ELL_PLUS, ELL_MINUS = (Q(1, 2), Q(1, 2)), (Q(1, 2), Q(-1, 2))


def odd_readouts(curvature):
    """Two oriented cross-sector channels, NOT an inner-product assumption."""
    curvature = matrix(curvature)
    if len(curvature) != 2 or mul(mul(K, curvature), K) != scale(curvature, -1):
        raise ValueError('This recovery contract requires a two-mode K-odd target')
    return (pair(ELL_PLUS, curvature, W_MINUS),
            pair(ELL_MINUS, curvature, W_PLUS))


def recover_odd(q_plus, q_minus):
    q_plus, q_minus = scalar(q_plus), scalar(q_minus)
    x = (q_plus - q_minus) / 2
    y = -(q_plus + q_minus) / 2
    return add(scale(R, x), scale(RK, y))


def transport_probe(covector, vector, change):
    change = matrix(change)
    back = inverse(change)
    if len(covector) != len(change) or len(vector) != len(change):
        raise ValueError('Probe and reference dimensions differ')
    left = tuple(sum((scalar(covector[i]) * back[i][j]
                      for i in range(len(change))), Q(0)) for j in range(len(change)))
    right = tuple(sum((change[i][j] * scalar(vector[j])
                       for j in range(len(change))), Q(0)) for i in range(len(change)))
    return left, right


@dataclass(frozen=True)
class ConnectionJet:
    """Values A_i and derivatives partials[i][j]=d_j A_i at one point.

    Directions are coordinate comparison directions; their count is distinct
    from carrier dimension. Gauge/reference changes here are constant only.
    A varying projector would require dP terms and is outside this adapter.
    """
    directions: tuple
    operators: tuple
    partials: tuple

    def __post_init__(self):
        names = tuple(self.directions)
        if not names or len(set(names)) != len(names) or any(not isinstance(x, str) or not x for x in names):
            raise ValueError('Declare distinct nonempty coordinate direction names')
        d = len(names)
        operators = same_carrier(*self.operators)
        if len(operators) != d or len(self.partials) != d or any(len(row) != d for row in self.partials):
            raise ValueError('A d-direction jet needs d operators and d by d partials')
        partials = tuple(tuple(matrix(a) for a in row) for row in self.partials)
        same_carrier(*operators, *(a for row in partials for a in row))
        object.__setattr__(self, 'directions', names)
        object.__setattr__(self, 'operators', operators)
        object.__setattr__(self, 'partials', partials)

    @classmethod
    def constant(cls, directions, operators):
        operators = same_carrier(*operators)
        z = scale(operators[0], 0)
        d = len(operators)
        return cls(tuple(directions), operators, tuple(tuple(z for _ in range(d)) for _ in range(d)))

    def curvature(self, i, j):
        d = len(self.directions)
        if not (0 <= i < d and 0 <= j < d):
            raise IndexError('Direction index outside the declared jet')
        return add(sub(self.partials[j][i], self.partials[i][j]),
                   commutator(self.operators[i], self.operators[j]))

    def components(self):
        return {(self.directions[i], self.directions[j]): self.curvature(i, j)
                for i in range(len(self.directions)) for j in range(i + 1, len(self.directions))}

    def pullback(self, jacobian, new_directions):
        """Constant affine base map x=J y; J has d rows and m columns."""
        names = tuple(new_directions)
        j = tuple(tuple(scalar(x) for x in row) for row in jacobian)
        d, m = len(self.directions), len(names)
        if not m or len(j) != d or any(len(row) != m for row in j):
            raise ValueError('The Jacobian must have d old rows and m new columns')
        z = scale(self.operators[0], 0)

        def total(terms):
            out = z
            for factor, a in terms:
                out = add(out, scale(a, factor))
            return out

        operators = tuple(total((j[i][a], self.operators[i]) for i in range(d)) for a in range(m))
        partials = tuple(tuple(total((j[i][a] * j[k][b], self.partials[i][k])
                                    for i in range(d) for k in range(d))
                               for b in range(m)) for a in range(m))
        return ConnectionJet(names, operators, partials)

    def change_reference(self, change):
        change = matrix(change)
        same_carrier(*self.operators, change)
        inverse(change)
        return ConnectionJet(self.directions,
                             tuple(conjugate(a, change) for a in self.operators),
                             tuple(tuple(conjugate(a, change) for a in row) for row in self.partials))

    def compress(self, projector):
        p = matrix(projector)
        same_carrier(*self.operators, p)
        if mul(p, p) != p:
            raise ValueError('A constant idempotent projector is required')
        compress = lambda a: mul(mul(p, a), p)
        return ConnectionJet(self.directions, tuple(map(compress, self.operators)),
                             tuple(tuple(map(compress, row)) for row in self.partials))


def linear_response_jet(response, forces):
    """Scalar-sector A_i=(LX)_i; F_ij=L_ji-L_ij for constant L."""
    l = matrix(response)
    x = tuple(map(scalar, forces))
    if len(x) != len(l):
        raise ValueError('Force vector dimension must match constant response matrix')
    values = tuple(matrix([[sum((l[i][j] * x[j] for j in range(len(x))), Q(0))]])
                   for i in range(len(x)))
    partials = tuple(tuple(matrix([[l[i][j]]]) for j in range(len(x))) for i in range(len(x)))
    return ConnectionJet(tuple(f'X{i}' for i in range(len(x))), values, partials)


def derivative_cancellation_jet(t):
    """At any rational (s,t), A=g^-1 dg for g=(I+sE12)(I+tE21)."""
    t = scalar(t)
    a = matrix([[t, 1], [-t*t, -t]])
    b = matrix([[0, 0], [1, 0]])
    dt_a = matrix([[1, 0], [-2*t, -1]])
    z = scale(a, 0)
    return ConnectionJet(('s', 't'), (a, b), ((z, dt_a), (z, z)))


def exact_witnesses():
    """Deterministic explicit examples; general proofs are in the R8 document."""
    p = matrix([[1, 0], [0, 0]])
    p_plus, p_minus = parity_projectors(K)
    f = commutator(K, R)
    a3, b3 = matrix([[1, 0, 0], [0, 0, 0], [0, 0, 0]]), matrix([[0, 0, 0], [0, 1, 0], [0, 0, 0]])
    p3 = sub(identity(3), scale(matrix([[1]*3 for _ in range(3)]), Q(1, 3)))
    directions = ConnectionJet.constant(('radial', 'tangent_1', 'tangent_2'), (K, R, RK))
    pure = derivative_cancellation_jet(Q(2, 3))
    # L(x,y)=diag(y,0) is symmetric, but omega=xy dx has d omega=-x dx^dy.
    x, y = Q(2), Q(3)
    symmetric_variable = ConnectionJet(('x', 'y'), (matrix([[x*y]]), matrix([[0]])),
                                      ((matrix([[y]]), matrix([[x]])), (matrix([[0]]), matrix([[0]]))))
    return {
        'KIR': {'curvature': f, 'trace': trace(f), 'active_cut_flips_sign': mul(mul(K, f), K) == scale(f, -1),
                'odd_readouts': odd_readouts(f), 'recovered': recover_odd(*odd_readouts(f)),
                'same_sector_plus': mul(mul(p_plus, f), p_plus),
                'same_sector_minus': mul(mul(p_minus, f), p_minus),
                'cross_sector_plus_minus': mul(mul(p_plus, f), p_minus)},
        'full_curved_reduced_flat': compression_report(K, R, p),
        'full_flat_reduced_curved': compression_report(a3, b3, p3),
        'three_direction_components': directions.components(),
        'one_direction_pullback_components': directions.pullback(((1,), (2,), (3,)), ('path',)).components(),
        'noncommuting_flat_connection': {'commutator': commutator(*pure.operators), 'curvature': pure.curvature(0, 1)},
        'constant_onsager': {'symmetric_curvature': linear_response_jet(((2, 1), (1, 3)), (4, 5)).curvature(0, 1),
                            'asymmetric_curvature': linear_response_jet(((2, 1), (4, 3)), (4, 5)).curvature(0, 1)},
        'state_dependent_symmetric_response': {'point': (x, y), 'response': matrix([[y, 0], [0, 0]]),
                                             'curvature': symmetric_variable.curvature(0, 1)},
    }
