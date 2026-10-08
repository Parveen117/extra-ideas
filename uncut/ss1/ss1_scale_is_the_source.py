"""SS1: the source of non-commutation is a scale. Exact sympy."""
import json, os, sympy as sp
S, V, b, a_, bw, gam = sp.symbols('S V b a b_w gamma', positive=True)
G = sp.Function('G')
out = {}
z = lambda e: sp.simplify(e.doit()) == 0
DS = lambda f: S*sp.diff(f, S); DV = lambda f: V*sp.diff(f, V); Db = lambda f: b*sp.diff(f, b)
# centre potential with one built-in scale b, covariant when state and scale are stretched together
U = S**(1/a_)*G(V/S**(bw/a_), b/S**(gam/a_))
out['T0_covariant'] = z(a_*DS(U) + bw*DV(U) + gam*Db(U) - U)
P = -sp.diff(U, V); eps = DS(P)/DV(P)
# T1 (the law): the combination that CP1 found to vanish is exactly the response of eps to the scale
out['T1_law'] = z(a_*DS(eps) + bw*DV(eps) + gam*Db(eps))
# T2: no dependence on the scale -> CP1's tie (one defect); and eps free of everything -> commuting (ON1-T2)
U0 = S**(1/a_)*G(V/S**(bw/a_), 1)
P0 = -sp.diff(U0, V); e0 = DS(P0)/DV(P0)
out['T2_no_scale_gives_CP1'] = z(a_*DS(e0) + bw*DV(e0))
# T3: ON1's witness P = S/(V - b): weights (0, 1, 1): D_V eps = - D_b eps, the defect IS the scale response
Pw = S/(V - b); ew = DS(Pw)/DV(Pw)
out['T3_witness'] = z(DV(ew) + Db(ew)) and sp.simplify(Db(ew)) != 0 and sp.limit(Db(ew), b, 0) == 0
# T4: the outer scale. XU3's defect depends on r_b and h only through x = h r_b: with l = 1/h the same law
rb, h = sp.symbols('r_b h', positive=True); f = sp.Function('f')
out['T4_outer_scale'] = z(rb*sp.diff(f(h*rb), rb) - h*sp.diff(f(h*rb), h))
# T5: two scales. Inner b and outer l: defects are functions of the pure ratios; with neither, none
l = sp.symbols('l', positive=True)
H = lambda p, q, r_: 1/(1 + p + q) + p**2*r_ + q*r_**2 + p*q
U2 = (S**(1/a_)*H(V/S**(bw/a_), b/S**(gam/a_), l/S**(gam/a_))).subs({a_: 2, bw: 3, gam: 1})
P2 = -sp.diff(U2, V); e2 = DS(P2)/DV(P2)
out['T5_two_scales'] = z(2*DS(e2) + 3*DV(e2) + (Db(e2) + l*sp.diff(e2, l))) and sp.simplify(Db(e2)) != 0 and sp.simplify(l*sp.diff(e2, l)) != 0
out = {k_: (bool(x_) if not isinstance(x_, str) else x_) for k_, x_ in out.items()}
out['pass'] = all(x_ for x_ in out.values() if isinstance(x_, bool))
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'SS1_RESULT.json'), 'w'), indent=1)
if __name__ == '__main__': print(out)
