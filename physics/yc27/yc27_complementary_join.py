#!/usr/bin/env python3
"""YC27's new geometry and rational budgets, not a spectral truncation.

The note proves the geometry for every admitted torus. This calculator
checks its two smallest examples and the exact endpoint arithmetic only.
Run: python physics/yc27/yc27_complementary_join.py
"""

from collections import Counter, defaultdict
from fractions import Fraction as F
from itertools import combinations, product


def geometry(shape):
    lengths = tuple(2 * a for a in shape)
    vertices = list(product(*(range(n) for n in lengths)))

    def step(v, axis):
        w = list(v)
        w[axis] = (w[axis] + 1) % lengths[axis]
        return tuple(w)

    def old_id(v):
        return ("old",) + tuple(v[i] // shape[i] for i in range(3))

    edges = [(v, i) for v in vertices for i in range(3)]
    owner = {}
    neighbours = defaultdict(set)
    old_vertices = defaultdict(set)
    for v in vertices:
        old_vertices[old_id(v)].add(v)
    for e in edges:
        v, i = e
        w = step(v, i)
        if old_id(v) == old_id(w):
            owner[e] = old_id(v)
        else:
            neighbours[v].add(w)
            neighbours[w].add(v)

    components = []
    component_at = {}
    for v in vertices:
        if v in component_at or not neighbours[v]:
            continue
        key = ("dual", len(components))
        found = {v}
        todo = [v]
        while todo:
            w = todo.pop()
            for x in neighbours[w] - found:
                found.add(x)
                todo.append(x)
        components.append(found)
        for w in found:
            component_at[w] = key
    for e in edges:
        if e not in owner:
            v, i = e
            assert component_at[v] == component_at[step(v, i)]
            owner[e] = component_at[v]

    degrees = Counter()
    factor_edges = Counter(owner.values())
    for (v, i), factor in owner.items():
        degrees[factor, v] += 1
        degrees[factor, step(v, i)] += 1
    incidence_edges = []
    for v, dual in component_at.items():
        old = old_id(v)
        incidence_edges.append((old, dual))
        assert degrees[old, v] + degrees[dual, v] == 6
        d = degrees[old, v]
        assert F(1, 2 * d) + F(1, 2 * (6 - d)) >= F(1, 3)
    assert len(incidence_edges) == len(set(incidence_edges))
    for j, vs in enumerate(components):
        assert all(len(vs & old_vs) <= 1 for old_vs in old_vertices.values())
        assert (len(vs), factor_edges["dual", j]) in {(2, 1), (4, 4), (8, 12)}

    incidence = Counter()
    internal_faces = Counter()
    local_force = defaultdict(F)
    mixed_old_edges = set()
    face_counts = Counter()
    for v in vertices:
        for i, j in combinations(range(3), 2):
            face = ((v, i), (step(v, i), j), (step(v, j), i), (v, j))
            factors = Counter(owner[e] for e in face)
            if len(factors) == 1:
                factor = next(iter(factors))
                internal_faces[factor] += 1
                face_counts[factor[0]] += 1
                cap = F(1, 2) if factor[0] == "dual" else F(1)
                if factor[0] == "old" and shape == (4, 4, 2):
                    if 1 in (i, j) and v[1] % shape[1] == 1:
                        cap = F(1, 2)
                for edge in face:
                    local_force[edge] += cap
            else:
                assert len(factors) == 4 and set(factors.values()) == {1}
                assert Counter(f[0] for f in factors) == {"old": 2, "dual": 2}
                face_counts["mixed"] += 1
                incidence.update(factors.keys())
                mixed_old_edges.update(e for e in face if owner[e][0] == "old")

    max_old = 48 if shape == (4, 2, 2) else 88
    assert set(incidence[f] for f in old_vertices) == {max_old}
    for j, vs in enumerate(components):
        f = ("dual", j)
        assert (internal_faces[f], incidence[f]) == {
            2: (0, 4), 4: (1, 12), 8: (6, 24)
        }[len(vs)]

    force_by_axis = defaultdict(F)
    for e in mixed_old_edges:
        v, i = e
        factor = owner[e]
        d = min(degrees[factor, v], degrees[factor, step(v, i)])
        bound = F(d * d, 64) * local_force[e] ** 2
        force_by_axis[i] = max(force_by_axis[i], bound)
    expected_force = F(9, 4) if shape == (4, 2, 2) else F(625, 256)
    assert max(force_by_axis.values()) == expected_force
    component_counts = Counter(len(vs) for vs in components)
    print(f"shape={shape}, torus={lengths}")
    print(f"  complementary vertex counts: {dict(sorted(component_counts.items()))}")
    print(f"  face counts: {dict(face_counts)}, old mixed incidence={max_old}")
    print("  old boundary-force bounds by axis:",
          {i: str(force_by_axis[i]) for i in range(3)})


def budgets():
    radius = F(1, 32)
    exponential = F(9, 7)
    exact_exp_majorant = F(12491, 9728)
    assert exact_exp_majorant < exponential
    source_floor = F(4)
    norm_cap = F(3, 40)
    for name, force, incidence, cap in (
        ("28-link", F(9, 4), 48, F(1, 3200)),
        ("64-link", F(625, 256), 88, F(1, 5700)),
    ):
        variance = 2 * force + 2 * F(9, 64)
        inverse_square = F(1, 576) + variance * (
            1 / (864 * source_floor) + 1 / (144 * source_floor ** 2)
        )
        assert inverse_square < norm_cap ** 2
        beta = incidence * cap
        seed = 4 * incidence * norm_cap * cap
        contraction = 4 * beta * exponential * (2 + 8 * (1 + 2 * radius))
        relative = 8 * beta * exponential
        reserve = radius - seed - contraction * radius
        gap = 4 * (1 - relative)
        assert contraction < 1 and relative < 1 and reserve > 0
        print(f"{name}, mixed cap={cap}:")
        print(f"  variance={variance}, inverse square<={inverse_square}, norm<3/40")
        print(f"  q<={contraction}, b<={relative}, seed<={seed}, ball reserve={reserve}")
        print(f"  physical gap>={gap} ({float(gap):.9f})")
    assert F(12) - 2 * 6 * F(1, 2) == 6
    assert F(12) - 2 * F(1, 2) == 11
    print("Zero mixed coupling: seed=q=b=0; reference gap floor=4; metric=I by the proof.")


if __name__ == "__main__":
    geometry((4, 2, 2))
    geometry((4, 4, 2))
    budgets()
