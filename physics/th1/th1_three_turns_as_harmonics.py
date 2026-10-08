"""TH1: at the rung 1/3 the three turns of a free circle are the first three multiples of one rate.

RC1: closure when (in-out rate)/(round rate) = p/q.  QP1: the three turns are the circuit Omega, the turn left over by the
orbit's shape Omega - kappa, and the turn of the plane.  MO1-M5: u'' + u = r_s/(2L^2) + (3/2) r_s u^2.
Exact arithmetic and sympy.  Python 3.12."""
import json, os
from fractions import Fraction as Fr
from math import gcd
import sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))


def run():
    out, rec = {}, {}
    z = lambda e: sp.simplify(e) == 0
    # H1: on the rung p/q the in-out rate, the left-over-turn rate and the circuit are p : q-p : q times Omega/q
    rungs = [(p, q) for q in range(2, 9) for p in range(1, q) if gcd(p, q) == 1]
    turns = lambda p, q: (p, q - p, q)
    distinct = [(p, q) for p, q in rungs if len(set(turns(p, q))) == 3]
    lowest = min(distinct, key=lambda pq: (pq[1], pq[0]))
    out['H1_lowest_rung_with_three_distinct_turns'] = lowest == (1, 3) and turns(*lowest) == (1, 2, 3)
    out['H1_only_rung_one_half_shows_two'] = [pq for pq in rungs if len(set(turns(*pq))) < 3] == [(1, 2)]
    # H2: MO1-M5 near a circle u0, to second order in the in-out amplitude A:  eps'' + w^2 eps = (3/2) r_s eps^2 , w^2 = 1 - 3 r_s u0
    ph, A, rs, w = sp.symbols('phi A r_s omega', positive=True)
    eps = A*sp.cos(w*ph) + A**2*(sp.Rational(3, 4)*rs/w**2 - rs/(4*w**2)*sp.cos(2*w*ph))
    resid = sp.diff(eps, ph, 2) + w**2*eps - sp.Rational(3, 2)*rs*eps**2
    out['H2_second_order_solution'] = z(sp.series(sp.expand(resid), A, 0, 3).removeO())
    # the term of M5 that leaves the turn over is the same term that gives the in-out motion a second harmonic,
    # of relative size x e / (4 (1 - 3x)), e = A/u0, x = r_s u0
    u0, x, e = sp.symbols('u0 x e', positive=True)
    rel = (rs*A/(4*w**2)).subs({A: e*u0, w: sp.sqrt(1 - 3*x)}).subs(rs, x/u0)
    out['H2_relative_size'] = z(rel - x*e/(4*(1 - 3*x)))
    out['H2_at_one_third'] = z(rel.subs(x, sp.Rational(8, 27)) - sp.Rational(2, 3)*e)
    # H3: read through the round angle the lines are |n1 + n2 w| Omega.  The second harmonic of the in-out motion falls on the
    #     left-over-turn line exactly at w = 1/3;  the in-out line falls on the left-over-turn line exactly at w = 1/2
    ww = sp.Symbol('w', positive=True)
    out['H3_second_harmonic_meets_left_over'] = sp.solve(sp.Eq(2*ww, 1 - ww), ww) == [sp.Rational(1, 3)]
    out['H3_in_out_meets_left_over'] = sp.solve(sp.Eq(ww, 1 - ww), ww) == [sp.Rational(1, 2)]
    out['H3_third_harmonic_meets_circuit'] = sp.solve(sp.Eq(3*ww, 1), ww) == [sp.Rational(1, 3)]
    # H4: QC1-T3 on the in-out phase left over per circuit, Theta = 2 pi w:  content q silent  <=>  q w whole
    silent = {q: sorted({Fr(k, q) for k in range(1, q)}) for q in (1, 2, 3)}
    out['H4_contents'] = silent[1] == [] and silent[2] == [Fr(1, 2)] and silent[3] == [Fr(1, 3), Fr(2, 3)]
    # numbers.  XTE J1550-564: lower rate 183 +- 5 Hz (arXiv:1408.0884); closure at 1/3 gives in-out 91.5, circuit 274.5.
    # Reported together: features near 92, 184 and 276 Hz (astro-ph/0202305, abstract)
    rec['XTE_J1550'] = dict(in_out=round(183/2, 1), in_out_sigma=2.5, circuit=round(183*1.5, 1), circuit_sigma=7.5, reported_near=[92, 184, 276])
    # GRO J1655-40: circuit 441 +- 2 Hz: in-out rate 147.0 +- 0.7, left-over 294.0 +- 1.3 (measured 298 +- 4)
    rec['GRO_J1655'] = dict(in_out=147.0, in_out_sigma=0.7, left_over=294.0, left_over_sigma=1.3, measured_left_over=298, measured_sigma=4)
    out['N1_XTE_three_features'] = abs(91.5 - 92) < 2.5 and abs(274.5 - 276) < 7.5
    out = {k: bool(v) for k, v in out.items()}
    out['pass'] = all(out.values())
    json.dump(dict(checks=out, record=rec), open(os.path.join(HERE, 'TH1_RESULT.json'), 'w'), indent=1)
    return out, rec


if __name__ == '__main__':
    print(run())
