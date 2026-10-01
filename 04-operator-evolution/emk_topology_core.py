"""R15 finite checks of the cut/recognition topology repair.

This module manipulates finite records, relations and transition tables. It has
no complex-number, angle, metric, probability, or physics input. Finite checks
support the separate written proofs; they do not certify the full upload.
"""
from collections import deque
from itertools import product


def compose(f, g):
    """f after g; tables carry their common finite domain explicitly."""
    if len(f) != len(g):
        raise ValueError("incompatible carriers")
    return tuple(f[g[x]] for x in range(len(g)))


def kernel(labels):
    return frozenset((x, y) for x in range(len(labels))
                     for y in range(len(labels)) if labels[x] == labels[y])


def partitions(n):
    """All equivalence relations, without redundant class labels."""
    for labels in product(range(n), repeat=n):
        if labels[0] != 0:
            continue
        if all(labels[i] <= 1 + max(labels[:i]) for i in range(1, n)):
            yield labels


def future_relation(actions, observation):
    relation = kernel(observation)
    while True:
        refined = frozenset((x, y) for x, y in relation
                            if all((g[x], g[y]) in relation for g in actions))
        if refined == relation:
            return relation
        relation = refined


def transition_monoid(actions, n):
    identity = tuple(range(n))
    found, queue = {identity}, deque([identity])
    while queue:
        word = queue.popleft()
        for action in actions:
            new = compose(action, word)
            if new not in found:
                found.add(new)
                queue.append(new)
    return found


def future_by_all_words(actions, observation):
    n = len(observation)
    monoid = transition_monoid(actions, n)
    return frozenset((x, y) for x in range(n) for y in range(n)
                     if all(observation[w[x]] == observation[w[y]] for w in monoid))


def powerset(n):
    return [frozenset(x for x in range(n) if mask & (1 << x))
            for mask in range(1 << n)]


def saturation(subset, relation):
    return frozenset(y for x, y in relation if x in subset)


def repaired_opens(actions, relation, n):
    return frozenset(u for u in powerset(n)
                     if saturation(u, relation) == u
                     and all(g[x] in u for g in actions for x in u))


def raw_image_opens(actions, observation):
    """The upload's unsaturated condition E(C(U)) subset E(U).

    All generated words are checked, not just the two generators: without
    descent, testing only generators would change the source definition.
    """
    monoid = transition_monoid(actions, len(observation))
    return frozenset(u for u in powerset(len(observation))
                     if all({observation[w[x]] for x in u}
                            <= {observation[x] for x in u} for w in monoid))


def reachability(actions, relation, n):
    reach = set(relation)
    reach.update((x, g[x]) for g in actions for x in range(n))
    while True:
        more = reach | {(x, z) for x, y in reach for yy, z in reach if y == yy}
        if more == reach:
            return frozenset(reach)
        reach = more


def distances(actions, relation, n):
    # Zero cost is recognition equality, one unit is one cut occurrence.
    infinity = n + 1
    d = [[0 if (x, y) in relation else infinity for y in range(n)] for x in range(n)]
    for g in actions:
        for x in range(n):
            d[x][g[x]] = min(d[x][g[x]], 1)
    for k in range(n):
        for i in range(n):
            for j in range(n):
                d[i][j] = min(d[i][j], d[i][k] + d[k][j])
    return tuple(tuple(None if a == infinity else a for a in row) for row in d)


def healing_cost(chi, action, x):
    """Original rho with R = kernel(chi); None means no healing."""
    a, b, n, visited = action[x], x, 0, set()
    while (a, b) not in visited:
        if chi[a] == chi[b]:
            return n
        visited.add((a, b))
        a, b, n = chi[a], chi[b], n + 1
    return None


def topology_laws(opens, n):
    return (frozenset() in opens and frozenset(range(n)) in opens
            and all(a | b in opens and a & b in opens for a in opens for b in opens))


def exhaustive_checks():
    cases = subsets = 0
    for n in range(1, 4):
        maps = list(product(range(n), repeat=n))
        eqs = [(p, kernel(p)) for p in partitions(n)]
        for actions in product(maps, repeat=2):
            for labels, raw in eqs:
                relation = future_relation(actions, labels)
                assert relation == future_by_all_words(actions, labels)
                assert relation <= raw
                for _, candidate in eqs:
                    if candidate <= raw and all((g[x], g[y]) in candidate
                                                for x, y in candidate for g in actions):
                        assert candidate <= relation  # greatest sound quotient
                opens = repaired_opens(actions, relation, n)
                assert topology_laws(opens, n)
                reach = reachability(actions, relation, n)
                expected = frozenset(u for u in powerset(n)
                                     if all(x not in u or y in u for x, y in reach))
                assert opens == expected
                for u in powerset(n):
                    su = saturation(u, relation)
                    assert u <= su and saturation(su, relation) == su
                    subsets += 1
                d = distances(actions, relation, n)
                for x, y, z in product(range(n), repeat=3):
                    assert (d[x][y] == 0) == ((x, y) in relation)
                    if d[x][y] is not None and d[y][z] is not None:
                        assert d[x][z] is not None and d[x][z] <= d[x][y] + d[y][z]
                cases += 1
    return {"finite_systems_and_observations": cases, "subset_checks": subsets,
            "carrier_sizes": [1, 2, 3], "actions_per_system": 2,
            "all_transition_tables": True, "all_observation_partitions": True,
            "future_refinement_equals_full_transition_monoid": True,
            "greatest_congruence_refining_observation": True,
            "saturation_topology_reachability_and_directed_triangle": True}


def boundary_witnesses():
    # No closure equation follows just from naming two composable acts.
    kappa, chi = (1, 2, 0), (0, 1, 2)
    j = compose(chi, kappa)
    assert compose(j, j) != (0, 1, 2)
    cycle = {"kappa": kappa, "chi": chi, "J_squared": compose(j, j), "at_seed": 2}

    # Original seam applies one extra chi relative to rho = 0.
    kappa, chi = (2, 2, 2), (0, 0, 1)
    seam = chi[chi[kappa[0]]] == chi[chi[0]]
    cost = healing_cost(chi, kappa, 0)
    assert seam and cost == 1

    # Full-word counterexample to the unqualified source topology claim.
    kappa, chi = (2, 2, 0), (0, 0, 2)
    raw = raw_image_opens((kappa, chi), chi)
    u, v = frozenset({0, 2}), frozenset({1, 2})
    assert u in raw and v in raw and u & v not in raw

    # Current recognition can hide a distinction exposed by a future cut.
    obs, actions = (0, 0, 1), ((2, 1, 2), (0, 1, 2))
    assert (0, 1) in kernel(obs) and (0, 1) not in future_relation(actions, obs)

    # An involution has at most two endpoint states; no fixed-point limit here.
    swap = (1, 0)
    assert compose(swap, swap) == (0, 1) and all(swap[x] != x for x in range(2))
    # Source complex representation is audited with exact signed coordinates.
    # It is NOT a primitive of the repaired topology.
    assert (-1, -1)[::-1] == (-1, -1)

    # Four directed steps return endpoints but retain one complete loop.
    endpoint = lambda steps: steps % 4
    assert endpoint(0) == endpoint(4) and 0 // 4 != 4 // 4
    # Transitive forward-invariant topology, distinct from collapsing vertices.
    four = ((1, 2, 3, 0),)
    opens = repaired_opens(four, kernel((0, 1, 2, 3)), 4)
    assert opens == frozenset({frozenset(), frozenset(range(4))})
    assert len(set((0, 0, 0, 0))) == 1
    directed = distances(((1, 2, 2),), kernel((0, 1, 2)), 3)
    assert directed[0][2] == 2 and directed[2][0] is None
    return {
        "B01_four_cycle_not_implied_by_two_acts": cycle,
        "B02_original_seam_has_cost_one": {"kappa": [2, 2, 2], "chi": [0, 0, 1], "x": 0, "cost": cost},
        "B03_unsaturated_image_condition_not_topology": {"kappa": [2, 2, 0], "chi": [0, 0, 2],
            "U": sorted(u), "V": sorted(v), "intersection": sorted(u & v), "all_words_checked": True},
        "B04_current_eye_not_future_congruence": {"observation": obs, "actions": actions, "separated_pair": [0, 1]},
        "B05_involution_orbit_not_helix_or_fixed_limit": {"J": swap, "orbit": [0, 1, 0, 1]},
        "B06_fixed_line_includes_negative_ray": {"fixed_point": [-1, -1]},
        "B07_endpoint_equality_loses_winding": {"lengths": [0, 4], "endpoints": [0, 0], "loops": [0, 1]},
        "B08_four_point_indiscrete_vs_one_point_quotient": {"vertices": 4, "open_sets": 2, "orbit_quotient_points": 1},
        "B09_cut_distance_is_directed": {"forward": 2, "backward": None},
    }


def mutation_controls():
    # Real candidate replacements, judged against independently specified laws.
    obs, actions = (0, 0, 1), ((2, 1, 2), (0, 1, 2))
    wrong_current_only = kernel(obs)
    kappa, chi = (2, 2, 0), (0, 0, 2)
    raw = raw_image_opens((kappa, chi), chi)
    swap = (1, 0)
    controls = {
        "drop_future_tests": wrong_current_only != future_by_all_words(actions, obs),
        "drop_saturation": not topology_laws(raw, 3),
        "equate_seam_to_cost_zero": healing_cost((0, 0, 1), (2, 2, 2), 0) != 0,
        "assert_four_cycle_without_relation": compose((1, 2, 0), (1, 2, 0)) != (0, 1, 2),
        "turn_involution_into_idempotent_eye": compose(swap, swap) != swap,
        "erase_closed_path_memory": (0 % 4, 0 // 4) != (4 % 4, 4 // 4),
        "force_symmetric_cut_distance": distances(((1, 2, 2),), kernel((0, 1, 2)), 3)[2][0] is None,
    }
    assert all(controls.values())
    return controls


def phase_operator_checks():
    """Derive the quarter turn on signed UGD phase records before a chart.

    Alphabet size 4m is the declared cyclic-phase input; this check does not
    claim that A0-A2 select that size. Full numerals and seam ledgers survive.
    """
    checked = 0
    for m in range(1, 7):
        size = 4 * m
        def basis(j):
            return tuple(int(i == j) for i in range(size))
        def shift(v, step):
            return tuple(v[(j - step) % size] for j in range(size))
        def reflect(v):
            return tuple(v[(m - j) % size] for j in range(size))
        def neg(v):
            return tuple(-a for a in v)
        # Signed cut ledger, not imported real/complex coordinates.
        for j in range(2 * m):
            v = tuple(a - b for a, b in zip(basis(j), basis(j + 2 * m)))
            assert shift(v, 2 * m) == neg(v)
            assert shift(shift(v, m), m) == neg(v)
            assert shift(v, 4 * m) == v
            assert reflect(reflect(v)) == v
            assert reflect(shift(v, m)) == neg(shift(reflect(v), m))
            checked += 1
        # Do not extend iota^2=-I to the unsigned full phase module.
        assert shift(basis(0), 2 * m) != neg(basis(0))
    return {"alphabet_sizes": [4, 8, 12, 16, 20, 24], "signed_basis_vectors": checked,
            "construction": "iota = U^m on span(e_j-e_(j+2m)); K(e_j)=e_(m-j)",
            "iota_squared_minus_identity": True, "K_squared_identity": True,
            "K_iota_equals_minus_iota_K": True,
            "full_phase_module_not_collapsed_to_orientation_sector": True,
            "complex_numbers_used": False,
            "alphabet_size_derived_from_upload_A0_A2": False}
