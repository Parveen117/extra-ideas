"""Exact finite-matrix examples for R2 Aghora returns and projections.

All arithmetic is rational. General exponential-flow claims are proved in
AGHORA_RETURN_R2.md. This module evaluates nilpotent exponentials exactly
and a rational Cayley family, without approximating a matrix exponential.
"""

from fractions import Fraction as Q


def scalar(value):
    if not isinstance(value, (int, Q)):
        raise TypeError('Use integers or Fraction values')
    return Q(value)


def matrix(rows):
    result = tuple(tuple(scalar(x) for x in row) for row in rows)
    if not result or any(len(row) != len(result) for row in result):
        raise ValueError('A nonempty square matrix is required')
    return result


def identity(n):
    return matrix([[int(i == j) for j in range(n)] for i in range(n)])


def scale(a, factor):
    factor = scalar(factor)
    return tuple(tuple(factor * x for x in row) for row in a)


def add(a, b):
    if len(a) != len(b):
        raise ValueError('Matrix dimensions differ')
    return tuple(tuple(x + y for x, y in zip(ar, br)) for ar, br in zip(a, b))


def sub(a, b):
    return add(a, scale(b, -1))


def mul(a, b):
    if len(a) != len(b):
        raise ValueError('Matrix dimensions differ')
    n = len(a)
    return tuple(tuple(sum((a[i][k] * b[k][j] for k in range(n)), Q(0))
                       for j in range(n)) for i in range(n))


def transpose(a):
    return tuple(zip(*a))


def inverse(a):
    n = len(a)
    unit = identity(n)
    work = [list(a[i]) + list(unit[i]) for i in range(n)]
    for column in range(n):
        pivot = next((r for r in range(column, n) if work[r][column]), None)
        if pivot is None:
            raise ValueError('Matrix is singular')
        work[column], work[pivot] = work[pivot], work[column]
        divisor = work[column][column]
        work[column] = [x / divisor for x in work[column]]
        for row in range(n):
            if row != column:
                factor = work[row][column]
                work[row] = [x - factor * y for x, y in zip(work[row], work[column])]
    return tuple(tuple(row[n:]) for row in work)


def conjugate(a, change):
    return mul(mul(change, a), inverse(change))


def validate_source_relations(aghora, generator):
    if len(aghora) != len(generator):
        raise ValueError('A and G must act on the same space')
    if mul(aghora, aghora) != identity(len(aghora)):
        raise ValueError('A squared must be identity')
    if mul(mul(aghora, generator), aghora) != scale(generator, -1):
        raise ValueError('The source relation A G A = -G is required')


def nilpotent_flow(generator, parameter):
    """Exactly exp(tG) when G squared is zero."""
    unit = identity(len(generator))
    if mul(generator, generator) != scale(unit, 0):
        raise ValueError('This exact exponential requires G squared = 0')
    return add(unit, scale(generator, parameter))


def cayley_flow(generator, parameter):
    """(I+sG)(I-sG)^-1; a rational reversible family, not exp(sG)."""
    unit = identity(len(generator))
    scaled = scale(generator, parameter)
    return mul(add(unit, scaled), inverse(sub(unit, scaled)))


def return_operator(aghora, generator, flow):
    validate_source_relations(aghora, generator)
    if mul(mul(aghora, flow), aghora) != inverse(flow):
        raise ValueError('The chosen flow must satisfy A U A = U^-1')
    return mul(aghora, flow)


def parity_projectors(return_map):
    unit = identity(len(return_map))
    if mul(return_map, return_map) != unit:
        raise ValueError('The return must be an involution')
    return scale(add(unit, return_map), Q(1, 2)), scale(sub(unit, return_map), Q(1, 2))


def compressed_return(return_map, bond):
    """Return C=BRB, B-C^2, and BR(I-B)RB, with algebraic checks."""
    unit = identity(len(return_map))
    if mul(return_map, return_map) != unit:
        raise ValueError('The return must be an involution')
    if mul(bond, bond) != bond:
        raise ValueError('The bond must be a projection')
    compressed = mul(mul(bond, return_map), bond)
    defect = sub(bond, mul(compressed, compressed))
    excursion = mul(mul(mul(mul(bond, return_map), sub(unit, bond)), return_map), bond)
    return compressed, defect, excursion


def projection_for_scalar(value):
    """A rank-one idempotent B with B diag(1,-1) B = value * B.

    This is generally an oblique projection, not an orthogonal one.
    """
    value = scalar(value)
    row = ((1 + value) / 2, (1 - value) / 2)
    return matrix((row, row))


def orthogonal_example(parameter=Q(1, 3)):
    """Chosen rational example: full return +/-1; compressed value 4/5."""
    aghora = matrix(((1, 0), (0, -1)))
    generator = matrix(((0, 1), (-1, 0)))
    flow = cayley_flow(generator, parameter)
    returned = return_operator(aghora, generator, flow)
    bond = matrix(((1, 0), (0, 0)))
    return aghora, generator, flow, returned, bond
