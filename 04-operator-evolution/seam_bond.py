"""R3: an explicit two-mode seam bond and its exact rational returns.

The added rule is range(B)=L and kernel(B)=A L. It is not a consequence
of the source ratio alone. Metric compatibility is a further condition.
All inputs are integers or Fraction values; no floating-point fits occur.
"""

from fractions import Fraction as Q

from aghora_return import (
    add, cayley_flow, compressed_return, identity, inverse, matrix, mul,
    return_operator, scalar, scale, sub, transpose,
)


AGHORA = matrix(((1, 0), (0, -1)))


def _vector2(values):
    result = tuple(scalar(x) for x in values)
    if len(result) != 2:
        raise ValueError('A two-component state is required')
    return result


def apply(operator, state):
    """Apply an exact two-by-two matrix to an exact two-component state."""
    operator, state = matrix(operator), _vector2(state)
    if len(operator) != 2:
        raise ValueError('A two-by-two operator is required')
    return tuple(sum((operator[i][j] * state[j] for j in range(2)), Q(0))
                 for i in range(2))


def seam_constraint(gamma, state):
    """The typed scalar map q_gamma: R^2 -> R, q(x,y)=x-gamma*y."""
    gamma, (x, y) = scalar(gamma), _vector2(state)
    return x - gamma * y


def projector_from_line_and_involution(aghora, line_vector):
    """Project onto span(u) along span(Au), when these lines are transverse.

    This construction is covariant under any simultaneous invertible
    change of coordinates of A and u. A singular seam pair is rejected.
    """
    aghora, u = matrix(aghora), _vector2(line_vector)
    if len(aghora) != 2 or mul(aghora, aghora) != identity(2):
        raise ValueError('A two-dimensional involution is required')
    v = apply(aghora, u)
    basis = matrix(((u[0], v[0]), (u[1], v[1])))
    try:
        basis_inverse = inverse(basis)
    except ValueError as error:
        raise ValueError('The seam and its Aghora image must be transverse') from error
    return mul(mul(basis, matrix(((1, 0), (0, 0)))), basis_inverse)


def seam_geometry(gamma):
    """Construct B, K=[A,B], and the compatible metric representative H.

    H is derived only after requiring A to be an isometry and B to be
    self-adjoint. Its positive overall scale remains arbitrary.
    """
    gamma = scalar(gamma)
    if not gamma:
        raise ValueError('gamma=0 is not a transverse seam')
    u = (gamma, Q(1))
    v = apply(AGHORA, u)
    bond = projector_from_line_and_involution(AGHORA, u)
    structure = sub(mul(AGHORA, bond), mul(bond, AGHORA))
    metric = matrix(((1 / gamma**2, 0), (0, 1)))
    basis = matrix(((u[0], v[0]), (u[1], v[1])))
    return {'gamma': gamma, 'aghora': AGHORA, 'bond': bond,
            'complex_structure': structure, 'metric': metric,
            'seam_vector': u, 'mirror_vector': v, 'basis': basis}


def compatible_generator(gamma, omega):
    """All real A-odd, H-skew generators in this two-mode realization."""
    return scale(seam_geometry(gamma)['complex_structure'], scalar(omega))


def metric_adjoint(operator, metric):
    """H^-1 M^T H; callers supply an invertible metric."""
    operator, metric = matrix(operator), matrix(metric)
    return mul(mul(inverse(metric), transpose(operator)), metric)


def weighted_norm2(state, metric):
    """v^T H v; a squared norm when H is symmetric positive definite."""
    state = _vector2(state)
    weighted = apply(metric, state)
    return sum((x * y for x, y in zip(state, weighted)), Q(0))


def rational_seam_return(gamma, omega, parameter):
    """Exact Cayley return A(I+sG)(I-sG)^-1, not A exp(sG).

    Here G=omega*K and z=omega*s. The compressed coefficient and defect
    are computed from the matrices; their closed forms are tested apart.
    """
    geometry = seam_geometry(gamma)
    omega, parameter = scalar(omega), scalar(parameter)
    generator = scale(geometry['complex_structure'], omega)
    flow = cayley_flow(generator, parameter)
    returned = return_operator(AGHORA, generator, flow)
    compressed, defect, excursion = compressed_return(returned, geometry['bond'])
    # A rank-one projection has trace one, so trace(BRB) is its coefficient.
    mu = sum((compressed[i][i] for i in range(2)), Q(0))
    defect_coefficient = sum((defect[i][i] for i in range(2)), Q(0))
    return {**geometry, 'omega': omega, 'parameter': parameter,
            'generator': generator, 'flow': flow, 'returned': returned,
            'compressed': compressed, 'defect': defect, 'excursion': excursion,
            'mu': mu, 'defect_coefficient': defect_coefficient}
