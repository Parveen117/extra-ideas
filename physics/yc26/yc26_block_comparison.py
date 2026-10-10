"""Exact YC26 hidden-sector comparison and zero-specialization arithmetic.

The scalar sectors stand for full harmonic spaces, not truncated YM states.
Only the new proof's positivity, source, norm and endpoint arithmetic is run.
"""

from itertools import product

import sympy as s


Q = s.Rational
VERTICES = tuple(product(range(4), range(2)))
EVEN_MASKS = tuple(m for m in range(256) if m.bit_count() % 2 == 0)
FACE_MASKS = tuple(
    (1 << i) | (1 << j)
    for i, u in enumerate(VERTICES)
    for j, v in enumerate(VERTICES)
    if i < j and sum(abs(a - b) for a, b in zip(u, v)) == 1
)
ODD_COST = tuple(Q(4) if x in (0, 3) else Q(15, 4) for x, z in VERTICES)


def diagonal(mask):
    return sum((ODD_COST[j] for j in range(8) if mask >> j & 1), s.S.Zero) if mask else s.Integer(10)


def reflected(mask, flip_x, flip_z):
    return sum(
        ((mask >> j) & 1)
        << VERTICES.index((3 - x if flip_x else x, 1 - z if flip_z else z))
        for j, (x, z) in enumerate(VERTICES)
    )


def comparison_data():
    remaining = set(EVEN_MASKS)
    orbits = []
    while remaining:
        root = min(remaining)
        orbit = sorted({reflected(root, a, b) for a, b in product((0, 1), repeat=2)})
        orbits.append(orbit)
        remaining.difference_update(orbit)
    index = {m: i for i, orbit in enumerate(orbits) for m in orbit}
    adjacency = s.zeros(len(orbits))
    for i, orbit in enumerate(orbits):
        for face in FACE_MASKS:
            adjacency[i, index[orbit[0] ^ face]] += 1
    d = s.diag(*(diagonal(orbit[0]) for orbit in orbits))
    source = s.Matrix([Q(1, 2) if orbit[0] in FACE_MASKS else 0 for orbit in orbits])
    multiplicity = s.diag(*(len(orbit) for orbit in orbits))
    return orbits, index, adjacency, d, source, multiplicity


def comparison(t, z, data):
    orbits, index, adjacency, d, source, multiplicity = data
    matrix = d - z * s.eye(len(orbits)) - t * adjacency
    rhs = s.ones(len(orbits), 1).row_join(t * source)
    solution = matrix.inv() * rhs
    y, x = solution[:, 0], solution[:, 1]
    # Lift back and verify the actual 128-row equations. Reflection reduction
    # is a solving convenience, not a restriction of the positivity claim.
    for mask in EVEN_MASKS:
        i = index[mask]
        residual_y = (diagonal(mask) - z) * y[i] - t * sum(y[index[mask ^ f]] for f in FACE_MASKS)
        residual_x = (diagonal(mask) - z) * x[i] - t * sum(x[index[mask ^ f]] for f in FACE_MASKS)
        if residual_y != 1 or y[i] <= 0:
            raise ArithmeticError(f"Positive full-sector witness failed at mask {mask}")
        if residual_x != (t / 2 if mask in FACE_MASKS else 0) or x[i] < 0:
            raise ArithmeticError(f"Full-sector source inverse failed at mask {mask}")
    return (t * source.T * multiplicity * x)[0], 1 + (x.T * multiplicity * x)[0]


def parent_bounds(z):
    a = 8 - z
    c = 4 * (1 / (Q(32, 3) - z) + 1 / (16 - z) + 2 / (8 - z))
    b = 4 * (1 / (Q(32, 3) - z) ** 2 + 1 / (16 - z) ** 2 + 2 / (8 - z) ** 2)
    alpha = a - c
    return alpha, 1 / alpha, 1 + (1 + b) / alpha**2


def main():
    data = comparison_data()
    print(f"Interface: {len(VERTICES)} bridges, {len(FACE_MASKS)} faces, {len(EVEN_MASKS)} sectors; {len(data[0])} reflection orbits.")
    for delta, z in ((Q(27, 5), Q(4)), (Q(11), Q(9, 2))):
        alpha, sigma, _ = parent_bounds(z)
        reserve = delta - z - sigma
        if alpha <= 0 or reserve <= 0:
            raise ArithmeticError("Parent physical gap reserve failed")
        print(f"28-link block: delta={delta}, z={z}, alpha={alpha}, reserve={reserve}")
    _, _, parent_metric = parent_bounds(Q(1, 2))
    if parent_metric > Q(26, 25):
        raise ArithmeticError("Parent metric bound failed")
    print(f"28-link metric through z=1/2: {parent_metric} <= 26/25")
    zero_sigma, zero_metric = comparison(s.S.Zero, s.S.Zero, data)
    if (zero_sigma, zero_metric) != (0, 1):
        raise ArithmeticError("Zero specialization failed")
    print("64-link zero specialization: return=0, graph metric=I")
    rows = (
        (Q(1), Q(0), Q(77, 100), Q(10, 7), None),
        (Q(1), Q(1, 2), Q(109, 100), Q(101, 50), None),
        (Q(1), Q(1), Q(113, 50), None, Q(4)),
        (Q(1), Q(11, 10), Q(3), None, Q(9, 2)),
        (Q(3, 4), Q(5, 2), Q(4, 5), None, Q(4)),
        (Q(1, 2), Q(7, 2), Q(6, 25), None, Q(4)),
        (Q(1, 2), Q(4), Q(1, 3), None, Q(9, 2)),
    )
    for t, z, return_cap, metric_cap, delta in rows:
        sigma, metric = comparison(t, z, data)
        if sigma > return_cap or (metric_cap is not None and metric > metric_cap):
            raise ArithmeticError("Outward response or metric bound failed")
        reserve = None if delta is None else delta - z - return_cap
        if reserve is not None and reserve <= 0:
            raise ArithmeticError("64-link physical gap reserve failed")
        print(f"t={t}, z={z}: hidden positive, return<={return_cap}, metric cap={metric_cap}, gap reserve={reserve}")
    print(f"Inherited YC25 lattice floors: {Q(4)*Q(209,245)}, {Q(9,2)*Q(1033,1225)}")
    print(f"Old 1/4000 lattice window: {Q(9,2)*Q(4087,4375)}")
    print(f"Consistent incidence/mu: rectangles={Q(48)/Q(15,8)}, half-strength transverse blocks={Q(88)/Q(9,5)}")


if __name__ == "__main__":
    main()
