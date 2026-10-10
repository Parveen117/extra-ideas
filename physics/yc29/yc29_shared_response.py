#!/usr/bin/env python3
"""YC29's new harmonic arithmetic and unwrapped local residue counts.

No finite-spectrum approximation, predecessor regression or generic test suite.
"""

from collections import Counter
from fractions import Fraction as F
from itertools import combinations, product


def move(v, axis, amount=1):
    w = list(v)
    w[axis] += amount
    return tuple(w)


def edges(face):
    v, i, j = face
    return {(v, i), (move(v, i), j), (move(v, j), i), (v, j)}


def around(edge):
    v, i = edge
    result = set()
    for j in range(3):
        if j != i:
            axes = tuple(sorted((i, j)))
            result.add((v, *axes))
            result.add((move(v, j, -1), *axes))
    assert len(result) == 4
    return result


def kind(face, shape):
    v, i, j = face
    return sum(v[k] % shape[k] == shape[k] - 1 for k in (i, j))


def boundary(n, size):
    return int(n % size in (0, size - 1))


def pure_degree(face, shape):
    v, i, j = face
    k = 3 - i - j
    b = lambda axis, n: boundary(n, shape[axis])
    face_kind = kind(face, shape)
    assert face_kind != 1
    if face_kind == 0:
        return (b(i, v[i]) + b(i, v[i] + 1)
                + b(j, v[j]) + b(j, v[j] + 1) + 4 * b(k, v[k]))
    return 12 - 4 * b(k, v[k])


def neighbors(face, wanted, shape):
    candidates = set().union(*(around(e) for e in edges(face)))
    return {q for q in candidates if wanted(kind(q, shape))}


def centre_signs(face, shape):
    """The new all-face cut C composed with YC28's represented J."""
    def omega(v):
        x, y, z = v
        return (x * y + x * z + y * z) % 2

    def old(v):
        return tuple(v[k] // shape[k] for k in range(3))

    c = j = 0
    for v, i in edges(face):
        w = move(v, i)
        c ^= sum(v[:i]) % 2
        if old(v) == old(w):
            j ^= (omega(v) + omega(w)) % 2
    assert c == 1
    assert (c ^ j) == (kind(face, shape) != 1)


def local_count(shape, expected):
    rows = Counter()
    columns = Counter()
    for v in product(*(range(a) for a in shape)):
        for i, j in combinations(range(3), 2):
            p = (v, i, j)
            centre_signs(p, shape)
            if kind(p, shape) != 1:
                continue
            adjacent = neighbors(p, lambda k: k != 1, shape)
            assert len(adjacent) <= 8
            columns[len(adjacent)] += 1
            total = 0
            for r in adjacent:
                assert len(edges(p) & edges(r)) == 1
                actual = neighbors(r, lambda k: k == 1, shape)
                assert all(len(edges(r) & edges(q)) == 1 for q in actual)
                degree = pure_degree(r, shape)
                assert degree == len(actual)
                assert degree <= (8 if kind(r, shape) == 0 else 12)
                total += degree
            rows[total] += 1
    assert dict(sorted(rows.items())) == expected
    budget = max(rows)
    c = F(1, 22464)
    print(f"shape={shape}, mixed faces per cell={sum(rows.values())}")
    print(f"  pure neighbors per mixed face: {dict(sorted(columns.items()))}")
    print(f"  exact M* M row sums: {dict(sorted(rows.items()))}")
    print(f"  Gamma_2 norm <= {budget} / (48 * 22464^2) = {budget * c*c / 48}")
    print(f"  neutral norm kernel L_2 <= {budget * c*c / 576}")
    print("  C flips every face; D=CJ flips pure faces and fixes mixed faces")


def harmonic_arithmetic():
    shared = (F(1, 18) + F(3, 26)) / 16
    disjoint = F(1, 96)
    c = 2 * shared / 12 - F(1, 576)
    assert shared == F(5, 468)
    assert c == F(1, 22464)
    assert 2 * disjoint / 12 - F(1, 576) == 0

    # Independent heat-kernel evaluation of the same first source derivative.
    heat = (F(1, 84 * 12) - F(1, 48 * 18) + F(1, 112 * 26)) / 4
    assert heat == c
    gram = c * c * F(1, 4) / 12
    norm = c * c * F(1, 4) / (12 * 12)
    assert gram == 12 * norm
    print(f"Shared-edge inverse coefficient={shared}; source derivative={c}")
    print("No shared edge: source derivative exactly zero")
    print(f"One common internal plaquette: neutral p^2 q^2 coefficient={-2 * gram}")


if __name__ == "__main__":
    harmonic_arithmetic()
    local_count((4, 2, 2), {48: 4, 56: 8, 58: 4, 62: 8})
    local_count((4, 4, 2), {22: 4, 40: 8, 46: 8, 54: 16, 56: 4, 60: 4})
