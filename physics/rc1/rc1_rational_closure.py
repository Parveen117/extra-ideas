"""RC1: rational = closure, irrational = a turn left over.  The owner's statement, in the line's own formulas.

Exact rational arithmetic (fractions) for the arithmetic; sympy for the identities with symbols.  Python 3.12."""
import json, os
from fractions import Fraction as Fr
from math import gcd, isqrt
import sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))


def is_square(fr):
    n, d = fr.numerator, fr.denominator
    return n >= 0 and isqrt(n)**2 == n and isqrt(d)**2 == d


def run():
    out, rec = {}, {}
    pairs = [(p, q) for q in range(1, 13) for p in range(0, q) if gcd(p, q) == 1]
    # R1: MO1-M5: (in-out rate / round rate)^2 = 1 - 3x, x = r_s/r.  The orbit closes iff that ratio is rational p/q:
    #     r / r_s = 3 q^2 / (q^2 - p^2) ; it closes after q circuits.
    ok = True
    for p, q in pairs:
        rr = Fr(3*q*q, q*q - p*p)
        ok &= (1 - 3/rr == Fr(p, q)**2)
    out['R1_orbit_ladder'] = ok
    # OA1-T3: carried direction: (turn against the co-turning basis per circuit / 2 pi)^2 = 1 - 3x/2
    ok = True
    for p, q in pairs:
        rr = Fr(3*q*q, 2*(q*q - p*p))
        ok &= (1 - Fr(3, 2)/rr == Fr(p, q)**2)
    out['R1_direction_ladder'] = ok
    # R2: the zero ends of the two ladders are MO1's two special circles: last stable circle 3 r_s, light circle (3/2) r_s
    out['R2_zero_ends'] = (Fr(3*1, 1 - 0) == 3) and (Fr(3, 2*(1 - 0)) == Fr(3, 2))
    # between a ladder value and the next the ratio is irrational or of higher denominator: sample of non-ladder radii
    out['R2_generic_radius_not_closed'] = not any(is_square(1 - 3/Fr(n, 1)) for n in (5, 6, 7, 8, 9, 10, 11, 13))
    # R3: both close at once  <=>  1 - 3x = (P/Q)^2 and 1 - 3x/2 = (C/Q)^2  <=>  P^2 + Q^2 = 2 C^2
    #     <=>  ((Q+P)/2, (Q-P)/2, C) is a right triangle with whole sides: an exact turn of PT1, u = (A + B R)/C
    found = []
    for Q in range(1, 200):
        for P in range(0, Q):
            if gcd(P, Q) != 1:
                continue
            x = (1 - Fr(P, Q)**2)/3
            if is_square(1 - Fr(3, 2)*x):
                C = isqrt((P*P + Q*Q)//2)
                A, B = (Q + P)//2, (Q - P)//2
                assert (P + Q) % 2 == 0 and A*A + B*B == C*C and P*P + Q*Q == 2*C*C
                found.append((P, C, Q, str(1/x)))
    out['R3_doubly_closed_are_right_triangles'] = len(found) > 10
    # and the other way: every right triangle with whole sides gives one
    ok, radii = True, {}
    for m in range(2, 12):
        for n in range(1, m):
            if gcd(m, n) != 1 or (m - n) % 2 == 0:
                continue
            A, B, C = m*m - n*n, 2*m*n, m*m + n*n
            Q, P = A + B, abs(A - B)
            x = (1 - Fr(P, Q)**2)/3
            ok &= (1 - 3*x == Fr(P, Q)**2) and (1 - Fr(3, 2)*x == Fr(C, Q)**2) and Fr(1, 3) > x > 0
            # the orbit ratio is tan(45 degrees - phi), phi the angle of the exact turn: (1 - tan)/(1 + tan), tan = B/A
            ok &= (Fr(A - B, A + B) == (1 - Fr(B, A))/(1 + Fr(B, A)))
            radii[C] = str(1/x)
    out['R3_every_exact_turn_gives_one'] = ok
    rec['doubly_closed_first'] = found[:8]
    rec['radius_over_r_s_for_prime_turns'] = {str(k): radii[k] for k in (5, 13, 17, 29) if k in radii}
    out['R3_prime_5'] = radii.get(5) == '49/16'
    # R5: a point displaced by iota a, in n dimensions, is read where its distance vanishes: a sphere S^(n-2) of radius a
    a = sp.Symbol('a', positive=True)
    ok = True
    for nd in (2, 3, 4, 5):
        X = sp.symbols('X1:%d' % nd, real=True)
        Z = sp.Symbol('Z', real=True)
        R2 = sum(v**2 for v in X) + (Z - sp.I*a)**2
        re, im = sp.re(sp.expand(R2)), sp.im(sp.expand(R2))
        ok &= sp.simplify(im + 2*a*Z) == 0 and sp.simplify(re.subs(Z, 0) - (sum(v**2 for v in X) - a**2)) == 0
    out['R5_point_read_as_sphere'] = bool(ok)
    # R6: the shells of constant real part of that distance are confocal ellipses (ellipsoids); the foci are the read points
    r, th = sp.symbols('r theta', positive=True)
    Xc, Zc = sp.sqrt(r**2 + a**2)*sp.sin(th), r*sp.cos(th)
    out['R6_ellipse'] = sp.simplify(Xc**2/(r**2 + a**2) + Zc**2/r**2 - 1) == 0
    out['R6_complex_distance'] = sp.simplify(sp.expand(Xc**2 + (Zc - sp.I*a)**2 - (r - sp.I*a*sp.cos(th))**2)) == 0
    d_plus2 = sp.expand((Xc - a)**2 + Zc**2); d_minus2 = sp.expand((Xc + a)**2 + Zc**2)
    out['R6_foci'] = (sp.simplify(d_plus2 - (sp.sqrt(r**2 + a**2) - a*sp.sin(th))**2) == 0
                      and sp.simplify(d_minus2 - (sp.sqrt(r**2 + a**2) + a*sp.sin(th))**2) == 0)
    # R7: first-order orbits (MO1-M5 without its last term) are ellipses; their velocities, lifted to the 3-sphere in four
    #     dimensions, run along great circles: the lifted points stay in one plane through the centre
    k, h, e, ph = sp.symbols('k h e phi', positive=True)
    vx, vy = -(k/h)*sp.sin(ph), (k/h)*(e + sp.cos(ph))
    p0 = (k/h)*sp.sqrt(1 - e**2)                                   # p0^2 = -2 x energy
    v2 = vx**2 + vy**2
    xi = [2*p0*vx/(v2 + p0**2), 2*p0*vy/(v2 + p0**2), (v2 - p0**2)/(v2 + p0**2)]
    out['R7_on_the_sphere'] = sp.simplify(sum(c_**2 for c_ in xi) - 1) == 0
    out['R7_great_circle'] = sp.simplify(e*xi[1] - sp.sqrt(1 - e**2)*xi[2]) == 0
    u = (1 + e*sp.cos(ph))*k/h**2
    out['R7_orbit_equation'] = sp.simplify(sp.diff(u, ph, 2) + u - k/h**2) == 0
    out = {k_: bool(v) for k_, v in out.items()}
    out['pass'] = all(out.values())
    json.dump(dict(checks=out, record=rec), open(os.path.join(HERE, 'RC1_RESULT.json'), 'w'), indent=1)
    return out, rec


if __name__ == '__main__':
    print(run())
