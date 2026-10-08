"""XU2: (i) when the expansion rate varies the count needs a third plane, whose partner is a volume;
(ii) XU1's one-sign claim tested against reading the energy with the same frame. Exact sympy + sampled signs."""
import json, os, sympy as sp
rb, h = sp.symbols('r_b h', positive=True)
R_ = sp.Rational
out = {}
z = lambda e: sp.simplify(e) == 0
Mb = rb*(1 - h**2*rb**2)/2
kap = (1 - 3*h**2*rb**2)/(2*rb); P = 2*sp.pi/kap
S = sp.pi*rb**2
curl = lambda tr, th: sp.diff(th, rb) - sp.diff(tr, h)
# (i) the (M, t)-area alone is not the count when h varies
c0 = sp.simplify(curl(P*sp.diff(Mb, rb), P*sp.diff(Mb, h)))
out['T1_Mt_area_not_closed'] = z(c0 - 8*sp.pi*h*rb**3*(3*h**2*rb**2 - 2)/(3*h**2*rb**2 - 1)**2)
# the missing plane: dS = P dM + P (h r_b^3) dh ; partner of h^2 is r_b^3/2 (a volume)
out['T2_third_plane'] = z(sp.diff(S, rb) - P*sp.diff(Mb, rb)) and z(sp.diff(S, h) - (P*sp.diff(Mb, h) + P*h*rb**3))
out['T2_partner_is_volume'] = z(-sp.diff(Mb, h**2 if False else h)/(2*h) - rb**3/2)
# (ii) forms of the count read by the best frame, unit c
c = sp.sqrt(1 - 3*(h*Mb)**R_(2, 3))
forms = {
 'A_dS_over_c': (sp.diff(S, rb)/c, sp.diff(S, h)/c),                                   # XU1
 'F_frame_energy': (sp.diff(S, rb)/c - P*Mb*sp.diff(c, rb)/c**2, sp.diff(S, h)/c - P*Mb*sp.diff(c, h)/c**2),
}
pts = [(R_(1, 2), R_(1, 5)), (1, R_(1, 3)), (2, R_(1, 5)), (3, R_(1, 10)), (R_(1, 10), 5), (1, R_(1, 2)), (1, R_(1, 20))]
signs = {}
for name, (a, b) in forms.items():
    f = sp.lambdify((rb, h), curl(a, b), 'mpmath')
    signs[name] = [int(sp.sign(sp.N(f(*p), 20))) for p in pts]
out['T3_signs_A'] = str(signs['A_dS_over_c']); out['T3_signs_F'] = str(signs['F_frame_energy'])
out['T3_A_one_sign'] = len(set(signs['A_dS_over_c'])) == 1
out['T3_F_one_sign_opposite'] = set(signs['F_frame_energy']) == {1} and set(signs['A_dS_over_c']) == {-1}   # the sense reverses
out['T3_never_closed'] = all(s != 0 for v in signs.values() for s in v)
# (iii) with the outer horizon: total flat count is a function of state and falls as the centre grows
rc = (-rb + sp.sqrt(4/h**2 - 3*rb**2))/2
out['T4_outer_root'] = z((1 - 2*Mb/rc - h**2*rc**2))
St = sp.pi*(rb**2 + rc**2)
out['T4_total_falls'] = all(sp.N(sp.diff(St, rb).subs({rb: a, h: b})) < 0 for a, b in pts)
out['T4_empty_space_is_largest'] = z(St.subs(rb, 0) - sp.pi/h**2)
out = {k_: (bool(x_) if not isinstance(x_, str) else x_) for k_, x_ in out.items()}
out['pass'] = all(x_ for x_ in out.values() if isinstance(x_, bool))
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'XU2_RESULT.json'), 'w'), indent=1)
if __name__ == '__main__': print(out)
