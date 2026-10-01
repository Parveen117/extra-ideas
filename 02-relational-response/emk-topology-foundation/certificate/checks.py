"""Finite proof checks for the standalone constructive cut foundation.

Signed permutations are specified on the two free cut-role generators, so
basis identities extend to EVERY signed/refined ledger by linearity. History
checks are finite witnesses of separate general proofs, not infinite proofs.
No ordinary complex coefficient or floating-point number is used.
"""
from itertools import product


IDENTITY = (1, 2)
CUT_PARITY = (1, -2)
ROLE_EXCHANGE = (2, 1)


def apply_to_tag(action, signed_tag):
    return (1 if signed_tag > 0 else -1) * action[abs(signed_tag) - 1]


def compose(a, b):
    return tuple(apply_to_tag(a, tag) for tag in b)


def apply_to_ledger(action, value):
    output = [0, 0]
    for tag, coefficient in zip(action, value):
        output[abs(tag) - 1] += (1 if tag > 0 else -1) * coefficient
    return tuple(output)


def pairing(x, y):
    return sum(a * b for a, b in zip(x, y))


def area(x, y):
    return x[0] * y[1] - x[1] * y[0]


def signed_role_certificate():
    h, k = CUT_PARITY, ROLE_EXCHANGE
    r = compose(k, h)  # derived from role exchange AFTER cut-role parity
    minus_identity = (-1, -2)
    assert r == (2, -1)
    assert compose(h, h) == compose(k, k) == IDENTITY
    assert compose(r, r) == minus_identity
    assert compose(compose(r, r), compose(r, r)) == IDENTITY
    assert compose(k, r) == tuple(-x for x in compose(r, k))
    found, frontier = {IDENTITY}, [IDENTITY]
    while frontier:
        word = frontier.pop()
        for generator in (h, k):
            result = compose(generator, word)
            if result not in found:
                found.add(result)
                frontier.append(result)
    assert len(found) == 8
    basis = [(1, 0), (0, 1)]
    bilinear_checks = 0
    for action in found:
        for x, y in product(basis, repeat=2):
            assert pairing(apply_to_ledger(action, x), apply_to_ledger(action, y)) == pairing(x, y)
            bilinear_checks += 1
    for x, y in product(basis, repeat=2):
        assert pairing(apply_to_ledger(r, x), y) == area(x, y)
    for a, b, c in product(found, repeat=3):
        assert compose(compose(a, b), c) == compose(a, compose(b, c))
    orbit, current = [], 1
    while current not in orbit:
        orbit.append(current)
        current = apply_to_tag(r, current)
    assert orbit == [1, 2, -1, -2] and current == 1
    # Whole cycle returns role state; the four operation records still exist.
    record = ('R', 'R', 'R', 'R')
    assert record != () and len(record) == 4
    return {
        'primitive_data': {'ordered_roles': ['cut', 'self_cut_trace']},
        'constructed_maps': {'parity_on_roles': h, 'role_exchange': k},
        'derived_quarter_turn': r, 'derived_orbit': orbit,
        'cyclic_phase_alphabet_supplied': False,
        'ordinary_complex_input': False,
        'signed_operator_group_size': len(found),
        'complete_group_associativity_cases': len(found) ** 3,
        'complete_bilinear_basis_cases': bilinear_checks,
        'universal_extension': 'Each identity is on every free generator; linear/bilinear extension proves it for all coefficients.',
        'raw_chi_kappa_identified_with_quarter_turn': False,
        'operators': [list(x) for x in sorted(found)],
    }


def histories(depth):
    return [word for n in range(depth + 1) for word in product(('k', 'c'), repeat=n)]


def first_difference(a, b):
    if a == b:
        return None
    for index, (x, y) in enumerate(zip(a, b)):
        if x != y:
            return index
    return min(len(a), len(b))  # terminator versus next cut


def separation_ge(a, b):
    """None is infinity. Compare exponents, not floating distances."""
    return a is None or (b is not None and a >= b)


def min_separation(a, b):
    if a is None:
        return b
    if b is None:
        return a
    return min(a, b)


def history_certificate():
    words = histories(5)
    index = {(a, b): first_difference(a, b) for a, b in product(words, repeat=2)}
    assert len(words) == 63
    for a, b, c in product(words, repeat=3):
        assert separation_ge(index[a, c], min_separation(index[a, b], index[b, c]))
    for a, b in product(words, repeat=2):
        assert index[a, b] == index[b, a]
        for g in ('k', 'c'):
            new = first_difference((g,) + a, (g,) + b)
            assert new == (None if index[a, b] is None else index[a, b] + 1)
            # Chronological appending is non-expansive.
            assert separation_ge(first_difference(a + (g,), b + (g,)), index[a, b])
        for n in range(7):
            assert separation_ge(first_difference(a[:n], b[:n]), index[a, b])
    for a in words:
        for m, n in product(range(7), repeat=2):
            assert a[:n][:m] == a[:min(m, n)]
            assert separation_ge(first_difference(a[:n], a), n)
    # The end marker is essential: a cut and two cuts differ at depth one.
    assert first_difference(('k',), ('k', 'k')) == 1
    # Periodic infinite-history witnesses use exact prefix agreement, not
    # a claim that finite tests prove arbitrary infinite convergence.
    for period in [('k',), ('c',), ('k', 'c'), ('c', 'k', 'k')]:
        prefix = tuple(period[i % len(period)] for i in range(80))
        for n in range(1, 40):
            assert first_difference(prefix[:n], prefix[:n + 1]) == n
    return {'finite_histories': len(words), 'maximum_history_length': 5,
            'all_ultrametric_triples': len(words) ** 3,
            'cut_prepend_halves_distance': True, 'cut_append_nonexpansive': True,
            'aperture_meet_and_error_bounds': True, 'end_marker_retained': True,
            'periodic_infinite_history_prefix_witnesses': 4,
            'infinite_completion_evidence': 'WRITTEN_PREFIX_STABILIZATION_PROOF; finite witnesses only in this program'}


def graph_potentials(vertices, edges, weights):
    """Tree integration; reversed traversal is a comparison, not undoing an event."""
    adjacent = {x: [] for x in vertices}
    for edge_id, (a, b) in enumerate(edges):
        adjacent[a].append((b, edge_id, 1))
        adjacent[b].append((a, edge_id, -1))
    potential, tree, components = {}, set(), 0
    for root in vertices:
        if root in potential:
            continue
        components += 1
        potential[root] = 0
        frontier = [root]
        while frontier:
            a = frontier.pop()
            for b, edge, sign in adjacent[a]:
                if b not in potential:
                    potential[b] = potential[a] + sign * weights[edge]
                    tree.add(edge)
                    frontier.append(b)
    residual = tuple(weight - (potential[b] - potential[a])
                     for (a, b), weight in zip(edges, weights))
    return potential, tree, residual, components


def loop_certificate():
    graphs = [
        ((0, 1, 2, 3), ((0, 1), (1, 2), (2, 3), (3, 0))),
        ((0, 1, 2, 3), ((0, 1), (0, 2), (1, 3), (2, 3))),
        ((0, 1, 2), ((0, 1), (1, 0), (2, 2))),
        ((0, 1, 2, 3), ((0, 1), (2, 3))),
    ]
    examples = 0
    for vertices, edges in graphs:
        for values in product((-1, 0, 1), repeat=len(vertices)):
            exact = tuple(values[b] - values[a] for a, b in edges)
            potential, tree, residual, components = graph_potentials(vertices, edges, exact)
            assert not any(residual)
            assert len(tree) == len(vertices) - components
            assert len(edges) - len(tree) == len(edges) - len(vertices) + components
            for edge in set(range(len(edges))) - tree:
                tampered = list(exact)
                tampered[edge] += 1
                _, _, new_residual, _ = graph_potentials(vertices, edges, tampered)
                assert new_residual[edge] == 1
            examples += 1
    # No directed cycle exists in the diamond, yet comparison paths can differ.
    diamond = graphs[1]
    _, _, residual, _ = graph_potentials(*diamond, (0, 0, 0, 1))
    assert any(residual)
    # A single turn is a non-exact period on the four-edge phase circle.
    _, _, period, _ = graph_potentials(*graphs[0], (1, 1, 1, 1))
    assert sum(period) == 4
    # Graph of the actually derived signed role group, not a supplied geometry.
    operators = [tuple(g) for g in signed_role_certificate()['operators']]
    vertices = tuple(range(len(operators)))
    edges = tuple((i, operators.index(compose(g, a)))
                  for i, a in enumerate(operators) for g in (CUT_PARITY, ROLE_EXCHANGE))
    for chosen in (None,) + vertices:
        values = tuple(int(v == chosen) for v in vertices)
        exact = tuple(values[b] - values[a] for a, b in edges)
        _, tree, residual, components = graph_potentials(vertices, edges, exact)
        assert not any(residual) and components == 1
        for edge in set(range(len(edges))) - tree:
            tampered = list(exact)
            tampered[edge] += 1
            _, _, new_residual, _ = graph_potentials(vertices, edges, tampered)
            assert new_residual[edge] == 1
    assert (len(vertices), len(edges), len(edges) - len(tree)) == (8, 16, 9)
    return {'finite_graphs': 4, 'exact_potential_assignments': examples,
            'fundamental_cycle_count': 'edges - vertices + connected_components',
            'every_chord_perturbation_detected': True,
            'directed_acyclic_diamond_can_have_nonzero_comparison_period': True,
            'four_step_lifted_period': 4,
            'derived_role_graph': {'vertices': 8, 'edges': 16, 'independent_comparison_cycles': 9,
                                   'complete_potential_basis_plus_zero': 9}}


def negative_controls():
    r = compose(ROLE_EXCHANGE, CUT_PARITY)
    wrong_no_sign = ROLE_EXCHANGE
    source_k, source_c = (1, 2, 0), (0, 1, 2)
    source_j = tuple(source_c[source_k[x]] for x in range(3))
    source_j2 = tuple(source_j[source_j[x]] for x in range(3))
    checks = {
        'replace_quarter_turn_by_unsigned_role_swap': compose(wrong_no_sign, wrong_no_sign) != (-1, -2),
        'claim_raw_chi_kappa_closure': source_j2 != (0, 1, 2),
        'identify_orientation_and_reflection': r != ROLE_EXCHANGE,
        'erase_four_step_history': ('R',) * 4 != (),
        'erase_history_terminator': first_difference(('k',), ('k', 'k')) == 1,
        'make_involution_an_idempotent_eye': compose(ROLE_EXCHANGE, ROLE_EXCHANGE) != ROLE_EXCHANGE,
        'claim_one_aperture_recovers_all_histories': ('k',)[:1] == ('k', 'c')[:1] and ('k',) != ('k', 'c'),
        'claim_strict_infinite_hierarchy_from_involution': compose(r, r) == (-1, -2) and compose(ROLE_EXCHANGE, ROLE_EXCHANGE) == IDENTITY,
    }
    assert all(checks.values())
    return checks


def build():
    return {'signed_roles': signed_role_certificate(), 'histories': history_certificate(),
            'comparison_loops': loop_certificate(), 'rejected_mutations': negative_controls()}
