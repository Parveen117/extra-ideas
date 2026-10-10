#!/usr/bin/env python3
"""Necessary YC28 parity/geometry and exact return arithmetic only.

No harmonic diagonalization, predecessor regression or certificate generation.
The note proves the all-volume claims; these are its smallest two tori.
"""

from collections import Counter, defaultdict
from fractions import Fraction as F
from itertools import combinations, product


def geometry(shape):
    lengths = tuple(2 * a for a in shape)
    vertices = list(product(*(range(n) for n in lengths)))
    vertex_id = {v: i for i, v in enumerate(vertices)}

    def step(v, axis):
        w = list(v)
        w[axis] = (w[axis] + 1) % lengths[axis]
        return tuple(w)

    def old(v):
        return tuple(v[i] // shape[i] for i in range(3))

    def weight(v):
        x, y, z = v
        return (x * y + x * z + y * z) % 2

    def edge_sign(v, i):
        w = step(v, i)
        return (weight(v) + weight(w)) % 2 if old(v) == old(w) else 0

    def face_vertices(v, i, j):
        return (v, step(v, i), step(step(v, i), j), step(v, j))

    def face_type(v, i, j):
        crossings = sum(v[k] % shape[k] == shape[k] - 1 for k in (i, j))
        return ("old", "mixed", "dual")[crossings]

    faces = []
    for v in vertices:
        for i, j in combinations(range(3), 2):
            sign = (edge_sign(v, i) + edge_sign(step(v, i), j)
                    + edge_sign(step(v, j), i) + edge_sign(v, j)) % 2
            kind = face_type(v, i, j)
            assert sign == (kind == "mixed")
            if kind == "mixed":
                faces.append((v, i, j))
    face_id = {f: i for i, f in enumerate(faces)}
    signatures = [sum(1 << vertex_id[v] for v in face_vertices(*f)) for f in faces]
    assert len(signatures) == len(set(signatures))
    for f in faces:
        assert sum(weight(v) for v in face_vertices(*f)) % 2 == 1

    # A four-distinct-face relation is exactly a collision of pair XORs.
    pairs = defaultdict(list)
    for p, q in combinations(range(len(faces)), 2):
        pairs[signatures[p] ^ signatures[q]].append((p, q))
    quartets = set()
    for same_signature in pairs.values():
        for first, second in combinations(same_signature, 2):
            four = frozenset(first + second)
            assert len(four) == 4
            quartets.add(four)

    tubes = {}
    for v in vertices:
        six = [(w, i, j) for i, j in combinations(range(3), 2)
               for w in (v, step(v, 3 - i - j))]
        mixed = frozenset(face_id[f] for f in six if f in face_id)
        if not mixed:
            continue
        assert len(mixed) == 4
        pure = [f for f in six if f not in face_id]
        kinds = {face_type(*f) for f in pure}
        assert len(pure) == 2 and len(kinds) == 1
        tubes[mixed] = next(iter(kinds))
    assert quartets == set(tubes)
    incidence = Counter(p for tube in tubes for p in tube)
    assert set(incidence.values()) == {2}
    assert 2 * len(tubes) == len(faces)
    print(f"shape={shape}, torus={lengths}: {len(faces)} mixed faces")
    print("  exact sign flip on every face; distinct source signatures")
    print(f"  {len(tubes)} four-face returns, all elementary tubes: {dict(Counter(tubes.values()))}")
    print("  exactly two tubes per mixed face")


def arithmetic():
    for name, variance, expected in (
        ("28-link", F(153, 32), F(179, 6144)),
        ("64-link", F(661, 128), F(2197, 73728)),
    ):
        susceptibility = F(1, 48) + variance / (144 * 4)
        assert susceptibility == expected
        print(f"{name} vacuum face susceptibility: [1/48, {susceptibility}]")
    repeated = F(1, 4) / (F(4) * F(16, 3) * F(4))
    assert repeated == F(3, 1024)
    pair_bounds = {
        shared: (F(2, 1) / F(16, 3) if shared else F(0)) / 64
                + F(4, 64 * (8 - 2 * shared))
        for shared in range(3)
    }
    assert pair_bounds == {0: F(1, 128), 1: F(25, 1536), 2: F(11, 512)}
    tube = F(4, 4**3) + F(2, 4**2 * 8)
    derivative = F(12, 4**4) + F(4, 4**3 * 8) + F(2, 4**2 * 8**2)
    assert (tube, derivative) == (F(5, 64), F(29, 512))
    free_source = F(1, 64) * (F(16, 12 * 18 * 24) + F(8, 12 * 24 * 24))
    assert free_source == F(11, 165888)
    print("Fourth-order return bounds:")
    print(f"  repeated face={repeated}; repeated pairs={pair_bounds}")
    print(f"  tube={tube}; spectator-energy derivative={derivative}")
    print(f"  free old-ended or dual-ended tube source={free_source} times W_left W_right")


if __name__ == "__main__":
    geometry((4, 2, 2))
    geometry((4, 4, 2))
    arithmetic()
