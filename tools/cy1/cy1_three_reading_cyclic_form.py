"""CY1: the cyclic form of three readings.

Source object: the "cyclic tensor"  T_AB = alpha_A beta_B + beta_A gamma_B + gamma_A alpha_B  of the catalog
draft, with alpha = grad ln I, beta = grad ln J, gamma = grad ln K, and the tensor built from it,
    D_AB = T_AB/(IJK) - (1/2) h_AB ln(IJK) + w^2 delta_AB .
Sources read before building: extra-ideas R1 (return multipliers and ratios of equal-endpoint gains do not
depend on the reference), Publications CID-1 (the split of a response depends on the metric, the obstruction
does not).

Y1  Relabelling.  A cyclic relabelling of (alpha, beta, gamma) leaves T unchanged; exchanging two of them
    gives the transpose.  So the symmetric part is the same for all six orders and the alternating part
    changes sign with the order.
Y2  Symmetric part.   T + T^t = sigma sigma^t - alpha alpha^t - beta beta^t - gamma gamma^t ,
    sigma = alpha + beta + gamma.   Trace:  2 tr_h T = h(sigma, sigma) - h(alpha,alpha) - h(beta,beta) - h(gamma,gamma).
Y3  Alternating part. T - T^t = (alpha - gamma) ^ (beta - gamma)   (u ^ v = u v^t - v u^t).
    It is unchanged when one covector is added to all three, and it is zero exactly when the two
    differences are parallel.
Y4  For readings I, J, K:  T - T^t = d ln(I/K) ^ d ln(J/K).  A common factor (a change of reference
    I, J, K -> cI, cJ, cK) leaves it unchanged; the symmetric part T + T^t changes by
    2 (sigma delta^t + delta sigma^t) + 6 delta delta^t ,  delta = d ln c,
    which is zero only for delta = 0 or delta = -(2/3) sigma  (c = (IJK)^(-2/3)).
    For monomial readings in x, y it is  det(exponent differences) dln x ^ dln y.
Y5  The draft's D_AB is symmetric only where d ln(I/K) ^ d ln(J/K) = 0.  Its displayed law
    "T_AB + T_BC + T_CA = 0" is not an equation between objects of one type; read on triples of
    index values it fails for general readings.

Exact rational arithmetic.  Python 3.12, standard library only.
"""
from fractions import Fraction as F
from itertools import permutations
import json
import os
import random

HERE = os.path.dirname(os.path.abspath(__file__))
DIM = 4


def outer(u, v):
    return [[a*b for b in v] for a in u]


def madd(*ms):
    return [[sum(m[i][j] for m in ms) for j in range(len(ms[0][0]))] for i in range(len(ms[0]))]


def mscale(m, c):
    return [[c*x for x in row] for row in m]


def transpose(m):
    return [list(r) for r in zip(*m)]


def cyclic(a, b, c):
    return madd(outer(a, b), outer(b, c), outer(c, a))


def wedge(u, v):
    return madd(outer(u, v), mscale(outer(v, u), -1))


def vadd(u, v, s=1):
    return [a + s*b for a, b in zip(u, v)]


def parallel(u, v):
    return all(u[i]*v[j] == u[j]*v[i] for i in range(len(u)) for j in range(len(u)))


def form(h, u, v):
    return sum(h[i][j]*u[i]*v[j] for i in range(len(u)) for j in range(len(u)))


def trace_with(h, m):
    return sum(h[i][j]*m[i][j] for i in range(len(m)) for j in range(len(m)))


def is_zero(m):
    return all(x == 0 for row in m for x in row)


def is_sym(m):
    return m == transpose(m)


def dlog_monomial(expo, point):
    """gradient of ln(x^a y^b) at a rational point: (a/x, b/y)"""
    return [F(e)/p for e, p in zip(expo, point)]


def run():
    out = {}
    rng = random.Random(9102026)
    vec = lambda: [F(rng.randint(-6, 6), rng.randint(1, 4)) for _ in range(DIM)]
    triples = [(vec(), vec(), vec()) for _ in range(40)]
    zero = [[F(0)]*DIM for _ in range(DIM)]

    # Y1
    ok_c = ok_t = True
    for a, b, c in triples:
        t = cyclic(a, b, c)
        ok_c &= cyclic(b, c, a) == t and cyclic(c, a, b) == t
        ok_t &= cyclic(b, a, c) == transpose(t) and cyclic(a, c, b) == transpose(t) and cyclic(c, b, a) == transpose(t)
    out['Y1_cyclic_relabelling_leaves_T_unchanged'] = ok_c
    out['Y1_exchange_gives_the_transpose'] = ok_t
    ok = True
    for a, b, c in triples[:10]:
        syms = {str(madd(cyclic(*p), transpose(cyclic(*p)))) for p in permutations((a, b, c))}
        ok &= len(syms) == 1
    out['Y1_symmetric_part_is_the_same_for_all_six_orders'] = ok

    # Y2
    h = [[F(2 if i == j else (1 if abs(i - j) == 1 else 0)) for j in range(DIM)] for i in range(DIM)]
    ok_s = ok_tr = True
    for a, b, c in triples:
        t = cyclic(a, b, c)
        s = vadd(vadd(a, b), c)
        rhs = madd(outer(s, s), mscale(outer(a, a), -1), mscale(outer(b, b), -1), mscale(outer(c, c), -1))
        ok_s &= madd(t, transpose(t)) == rhs
        ok_tr &= 2*trace_with(h, t) == form(h, s, s) - form(h, a, a) - form(h, b, b) - form(h, c, c)
        ok_tr &= trace_with(h, t) == form(h, a, b) + form(h, b, c) + form(h, c, a)
    out['Y2_symmetric_part'] = ok_s
    out['Y2_trace'] = ok_tr

    # Y3
    ok_a = ok_shift = True
    for a, b, c in triples:
        t = cyclic(a, b, c)
        alt = madd(t, mscale(transpose(t), -1))
        ok_a &= alt == wedge(vadd(a, c, -1), vadd(b, c, -1))
        d = vec()
        t2 = cyclic(vadd(a, d), vadd(b, d), vadd(c, d))
        ok_shift &= madd(t2, mscale(transpose(t2), -1)) == alt
    out['Y3_alternating_part_is_the_wedge_of_the_differences'] = ok_a
    out['Y3_unchanged_by_a_common_shift'] = ok_shift
    a, b, c = triples[0]
    d = vadd(b, a, -1)
    col = cyclic(a, vadd(a, d), vadd(a, [F(-3, 2)*x for x in d]))            # three collinear readings
    out['Y3_zero_when_the_differences_are_parallel'] = is_sym(col)
    out['Y3_non_zero_otherwise'] = all(
        is_sym(cyclic(a, b, c)) == parallel(vadd(a, c, -1), vadd(b, c, -1)) for a, b, c in triples)

    # Y4: monomial readings in two variables at rational points
    ok_det = ok_ref = ok_sym_formula = ok_sym_zero_iff = True
    pts = [(F(2), F(3)), (F(1, 2), F(5)), (F(7, 3), F(1, 4))]
    for _ in range(30):
        ei, ej, ek = [(rng.randint(-4, 4), rng.randint(-4, 4)) for _ in range(3)]
        ec = (rng.randint(1, 3), rng.randint(-3, -1))                       # exponents of the common factor c
        det = (ei[0] - ek[0])*(ej[1] - ek[1]) - (ei[1] - ek[1])*(ej[0] - ek[0])
        for p in pts:
            al, be, ga = (dlog_monomial(e, p) for e in (ei, ej, ek))
            t = cyclic(al, be, ga)
            alt = madd(t, mscale(transpose(t), -1))
            ok_det &= alt == [[F(0), F(det)/(p[0]*p[1])], [-F(det)/(p[0]*p[1]), F(0)]]
            sh = lambda e: (e[0] + ec[0], e[1] + ec[1])
            t2 = cyclic(*(dlog_monomial(sh(e), p) for e in (ei, ej, ek)))
            ok_ref &= madd(t2, mscale(transpose(t2), -1)) == alt
            dl, sg = dlog_monomial(ec, p), vadd(vadd(al, be), ga)
            change = madd(mscale(madd(outer(sg, dl), outer(dl, sg)), 2), mscale(outer(dl, dl), 6))
            ok_sym_formula &= madd(t2, transpose(t2)) == madd(t, transpose(t), change)
            special = all(x == 0 for x in dl) or all(3*d == -2*s for d, s in zip(dl, sg))
            ok_sym_zero_iff &= is_zero(change) == special
    out['Y4_alternating_part_is_the_determinant_of_exponent_differences'] = ok_det
    out['Y4_unchanged_by_a_common_factor_of_the_three_readings'] = ok_ref
    out['Y4_change_of_the_symmetric_part'] = ok_sym_formula
    out['Y4_it_vanishes_only_for_two_special_factors'] = ok_sym_zero_iff
    # the special factor exists: I = x^-3, J = y^3, K = 1, c = x^2 y^-2 = (IJK)^(-2/3)
    p = pts[0]
    base = [dlog_monomial(e, p) for e in ((-3, 0), (0, 3), (0, 0))]
    moved = [dlog_monomial(e, p) for e in ((-1, -2), (2, 1), (2, -2))]
    tb, tm = cyclic(*base), cyclic(*moved)
    out['Y4_special_factor_leaves_both_parts_unchanged'] = (
        madd(tb, transpose(tb)) == madd(tm, transpose(tm)) and not is_sym(tb)
        and madd(tb, mscale(transpose(tb), -1)) == madd(tm, mscale(transpose(tm), -1)))

    # Y5
    p = pts[0]
    al, be, ga = dlog_monomial((1, 0), p), dlog_monomial((0, 1), p), dlog_monomial((0, 0), p)   # I = x, J = y, K = 1
    t = cyclic(al, be, ga)
    out['Y5_draft_tensor_is_not_symmetric_in_general'] = not is_sym(t)
    al2, be2, ga2 = dlog_monomial((1, 2), p), dlog_monomial((3, 6), p), dlog_monomial((-1, -2), p)  # one ratio only
    out['Y5_symmetric_exactly_on_dependent_ratios'] = is_sym(cyclic(al2, be2, ga2))
    # the displayed law read on triples of index values: for all i, j, k it would force T = 0 in every example
    ok = True
    for tr in triples[:10]:
        t = cyclic(*tr)
        holds = all(t[i][j] + t[j][k] + t[k][i] == 0 for i in range(DIM) for j in range(DIM) for k in range(DIM))
        ok &= holds == is_zero(t) and not holds
    out['Y5_displayed_cyclic_law_read_on_index_triples_fails'] = ok

    out['pass'] = all(v for v in out.values() if isinstance(v, bool))
    return out


if __name__ == '__main__':
    out = run()
    with open(os.path.join(HERE, 'CY1_RESULT.json'), 'w') as fh:
        json.dump(out, fh, indent=1, sort_keys=True)
        fh.write('\n')
    for k, v in out.items():
        print(f'{k}: {v}')
