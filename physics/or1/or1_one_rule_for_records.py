"""OR1: one rule for records.  Least record memory (GE2-T2) among the readings that share a record:
neighbours in repetition give closure and whole numbers; neighbours among histories give the stationary count.

Carrier: the rational block z = a + bR, R^2 = -1, det z = a^2 + b^2 (PT1); a turn is det = 1.
GE2-T2 on turns: record mean S = sum p_s u_s , memory M = 1 - det S = sum p_s det(u_s - S).
Exact rational arithmetic (fractions); sympy for the identities with symbols.  Python 3.12."""
import json, os
from fractions import Fraction as Fr
import sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))


class Z:
    """a + bR over the rationals."""
    def __init__(s, a, b=0): s.a, s.b = Fr(a), Fr(b)
    def __add__(s, o): return Z(s.a + o.a, s.b + o.b)
    def __sub__(s, o): return Z(s.a - o.a, s.b - o.b)
    def __mul__(s, o):
        if not isinstance(o, Z): return Z(s.a*o, s.b*o)
        return Z(s.a*o.a - s.b*o.b, s.a*o.b + s.b*o.a)
    def det(s): return s.a*s.a + s.b*s.b
    def __eq__(s, o): return s.a == o.a and s.b == o.b
    def __pow__(s, n):
        out = Z(1)
        for _ in range(n): out = out*s
        return out


ONE = Z(1)
U5, U13, QT = Z(Fr(3, 5), Fr(4, 5)), Z(Fr(5, 13), Fr(12, 13)), Z(0, 1)      # PT1's prime turns and the quarter turn


def mean(us, ps):
    S = Z(0)
    for u, p in zip(us, ps): S = S + u*p
    return S


def run():
    out, rec = {}, {}
    # A: GE2-T2 on turns: memory = 1 - det S = weighted variance; two readings: p(1-p) det(u_a - u_b)
    ok = True
    for us, ps in (([U5, U13, QT], [Fr(1, 2), Fr(1, 3), Fr(1, 6)]), ([ONE, U5, U5*U5, U13], [Fr(1, 4)]*4), ([U5, U13], [Fr(2, 7), Fr(5, 7)])):
        S = mean(us, ps)
        ok &= (1 - S.det() == sum(p*(u - S).det() for u, p in zip(us, ps)))
    out['A1_memory_is_variance'] = ok
    ok = True
    for ua, ub in ((U5, U13), (ONE, U5), (QT, U13)):
        for p in (Fr(1, 2), Fr(1, 3), Fr(3, 8)):
            ok &= (1 - mean([ua, ub], [1 - p, p]).det() == p*(1 - p)*(ua - ub).det())
    out['A2_two_readings_share_law'] = ok
    # B: repetition.  K successive returns of a turn u, equal weights.  Exact for every K:
    #    det( sum_{k<K} u^k ) . det(u - 1) = det(u^K - 1) <= 4
    ok, worst = True, Fr(0)
    for u in (U5, U13, U5*U13):
        acc, power = Z(0), ONE
        d1 = (u - ONE).det()
        for K in range(1, 301):
            acc = acc + power
            power = power*u
            ok &= (acc.det()*d1 == (power - ONE).det()) and (power - ONE).det() <= 4
            ok &= acc.det() <= 4/d1                                 # the sum never grows: bounded for every K
            worst = max(worst, acc.det()*d1)
        ok &= not any((u**n) == ONE for n in range(1, 60))         # PT1-P4: it never returns
    out['B1_turn_that_never_closes_has_bounded_sum'] = ok
    rec['largest_det_sum_times_det_u_minus_1'] = float(worst)
    # so the seen part of the record, det S_K, is at most 4 / (K^2 det(u-1)) for every K: explicit bound, memory -> 1
    out['B1_explicit_bound_u5'] = (4/((U5 - ONE).det()) == 5)
    # closure: the quarter turn in content 4 (u^4 = 1): the sum is K, S = 1, memory 0; in content 1 it is 0 after 4 returns
    acc1, acc4 = Z(0), Z(0)
    for k in range(8):
        acc1 = acc1 + QT**k
        acc4 = acc4 + (QT**4)**k
    out['B2_closed_content_grows'] = acc4 == Z(8) and acc1 == Z(0)
    # the cost between two successive returns is (1/4) det(u - 1) = (2 - u - conj u)/4 : CL1's clock curvature over 4
    for u in (U5, U13):
        out['B3_pair_cost_is_clock_curvature_%s' % ('5' if u is U5 else '13')] = Fr(1, 4)*(u - ONE).det() == (2 - 2*u.a)/4
    Th = sp.Symbol('Theta', real=True)
    out['B3_sin_squared'] = sp.simplify((2 - 2*sp.cos(Th))/4 - sp.sin(Th/2)**2) == 0
    # least: the pair cost is zero exactly at a whole number of turns, and those are its minima
    out['B4_least_at_closure'] = sp.solveset(sp.Eq(sp.sin(Th/2)**2, 0), Th, sp.Interval(0, 2*sp.pi)) == sp.FiniteSet(0, 2*sp.pi)
    # content q (QC1-T3): the same with u^q: a rung p/q of RC1 is closed for contents that are multiples of q and for no other
    def mean_content(p, q, c):
        return sp.simplify(sp.expand_complex(sum(sp.exp(2*sp.pi*sp.I*sp.Rational(p*c*k, q)) for k in range(q))/q))
    out['B5_rungs'] = (mean_content(1, 3, 1) == 0 and mean_content(1, 3, 2) == 0 and mean_content(1, 3, 3) == 1
                       and mean_content(1, 2, 1) == 0 and mean_content(1, 2, 2) == 1 and mean_content(2, 5, 5) == 1 and mean_content(2, 5, 3) == 0)
    # C: neighbours among histories.  Count tau(eps) = tau0 + a eps + b eps^2 ; reading Exp(g tau R) (CL2-T1: count of a leg = g tau)
    g, a, b, eps, t0 = sp.symbols('g a b epsilon tau0', real=True)
    tau = lambda e: t0 + a*e + b*e**2
    pair = sp.sin(g*(tau(eps) - tau(-eps))/2)**2                    # two-reading memory at equal share, x 1 (A2 with p = 1/2: det(u_a-u_b)/4)
    out['C1_pair_memory'] = sp.simplify(sp.series(pair, eps, 0, 3).removeO() - g**2*a**2*eps**2) == 0
    S3 = (sp.exp(sp.I*g*tau(-eps)) + sp.exp(sp.I*g*tau(0)) + sp.exp(sp.I*g*tau(eps)))/3
    M3 = sp.simplify(1 - sp.expand(S3*sp.conjugate(S3)))
    out['C2_three_neighbours'] = sp.simplify(sp.series(M3, eps, 0, 3).removeO() - sp.Rational(2, 3)*g**2*a**2*eps**2) == 0
    # so the memory among neighbouring histories vanishes to this order exactly when a = 0: the count is stationary
    # D: one function on MO1's circles (r_s = 1):  T = 2 pi sqrt(2) r^(3/2) , tau = T sqrt(1 - 3x/2) , E, L of MO1
    r = sp.Symbol('r', positive=True)
    x = 1/r
    T = 2*sp.pi*sp.sqrt(2)*r**sp.Rational(3, 2)
    tau_c = T*sp.sqrt(1 - sp.Rational(3, 2)*x)
    E = (1 - x)/sp.sqrt(1 - sp.Rational(3, 2)*x)
    L = sp.sqrt(r/(2 - 3*x))
    z = lambda e: sp.simplify(e) == 0
    out['D1_energy_is_slope_of_count'] = z(sp.diff(tau_c, r)/sp.diff(T, r) - E)
    W = E*T - tau_c
    out['D2_transform_is_the_area'] = z(W**2 - (2*sp.pi*L)**2) and W.subs(r, 5) > 0      # both sides positive outside the light circle
    out['D2_period_is_slope_of_area'] = z(sp.diff(W, r)/sp.diff(E, r) - T)
    # three closure quantities per circuit: clock phase g W ; carried direction 1 - tau/T (OA1) ; and their relation
    out['D3_turn_is_area_plus_binding'] = z((1 - tau_c/T) - (W/T + (1 - E)))
    xs = sp.Symbol('x', positive=True)
    first = sp.series(((W/T).subs(r, 1/xs)), xs, 0, 2).removeO(), sp.series((1 - E).subs(r, 1/xs), xs, 0, 2).removeO()
    out['D3_first_order_half_and_quarter'] = z(first[0] - xs/2) and z(first[1] - xs/4)
    out['D3_no_binding_on_the_diagonal'] = z(E.subs(r, 2) - 1)       # DG1: r = 2 r_s
    out = {k: bool(v) for k, v in out.items()}
    out['pass'] = all(out.values())
    json.dump(dict(checks=out, record=rec), open(os.path.join(HERE, 'OR1_RESULT.json'), 'w'), indent=1)
    return out, rec


if __name__ == '__main__':
    print(run())
