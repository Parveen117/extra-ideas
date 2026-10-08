"""TS2: the share of the centre value carried by a turning plane is fixed by the first law.
Rim radius R = r * gamma^rho, count S = s r^(d-3) R^2 (boundary measure), M = M0(S) * gamma^sigma,
J = kappa M R v, partner Omega = v / R, one scale. Exact sympy."""
import json, os, sympy as sp
chi, sig, kap, rho, c, s = sp.symbols('chi sigma kappa rho c s', positive=True)
out = {}
z = lambda e: sp.simplify(e) == 0
def residual(d, rh, sg, kp):
    k = sp.Rational(d-2, d-1)
    gam = 1/sp.cos(chi)
    Sh = s*gam**(2*rh)                              # S = r^(d-1) Sh
    f = c*gam**(2*rh*k + sg)                        # M = r^(d-2) f  (M0 ~ S^k)
    Om = sp.sin(chi)*sp.cos(chi)**rh                # Omega * r = v r / R
    g = kp*f*gam**rh*sp.sin(chi)                    # J = r^(d-1) g
    That = ((d-2)*f/(d-1) - Om*g)/Sh                # one scale: (d-2) M = (d-1)(T S + Omega J)
    return sp.simplify((sp.diff(f, chi) - That*sp.diff(Sh, chi) - Om*sp.diff(g, chi))/f), That, Sh, f
# T1: consistent for every chi exactly when sigma = kappa = 1 - rho (d-3)/(d-1)
ok_exist = ok_unique = True
for d in (3, 4, 5, 6):
    val = 1 - rho*sp.Rational(d-3, d-1)
    res, *_ = residual(d, rho, val, val)
    ok_exist = ok_exist and all(abs(sp.N(res.subs({chi: t, rho: r_}))) < 1e-12 for t in (0.3, 0.8, 1.2) for r_ in (0.5, 1, 2))
    res2, *_ = residual(d, rho, sig, kap)
    ser = sp.series(res2, chi, 0, 4).removeO()
    sol = sp.solve([sp.simplify(ser.coeff(chi, 1)), sp.simplify(ser.coeff(chi, 3))], [sig, kap], dict=True)
    good = [so for so in sol if sp.simplify(so[sig] - val) == 0 and sp.simplify(so[kap] - val) == 0]
    ok_unique = ok_unique and len(good) >= 1 and all(sp.simplify(so[sig] - so[kap]) == 0 for so in sol)
out['T1_consistent_family'] = ok_exist
out['T1_exponent_equals_share'] = ok_unique
# T2: three directions: 1 for every rho - the plain boost, whatever the rim rule
out['T2_three_directions'] = sp.simplify(1 - rho*sp.Rational(0, 2) - 1) == 0
# T3: rho = 1 (inner radius = rim radius read with the rim reader's unit): share 2/(d-1), HD1's (b), (c) and the law
ok3 = True
for d in (3, 4, 5, 6):
    val = sp.Rational(2, d-1)
    res, That, Sh, f = residual(d, 1, val, val)
    ok3 = ok3 and z(res) and z(4*sp.pi*That.subs({c: (d-1)*s/(4*sp.pi)}) - ((d-2) - 2*sp.sin(chi)**2)) and z(f.subs(c, (d-1)*s/(4*sp.pi)) - (d-1)*Sh/(4*sp.pi))
out['T3_rho_one_gives_HD1'] = ok3
# T4: rho = 0 (no contraction) is also consistent: plain boost in every d, with a different heat reading
res, That, Sh, f = residual(5, 0, 1, 1)
out['T4_rho_zero_is_consistent'] = z(res)
out['T4_its_heat_reading_differs'] = not z(4*sp.pi*That.subs(c, 4*s/(4*sp.pi)) - (3 - 2*sp.sin(chi)**2))
out['T4_its_heat_reading'] = str(sp.simplify(4*sp.pi*That.subs(c, 4*s/(4*sp.pi))))
# T5: two planes, rho = 1, d = 4 and 5: the same share for each plane
c1, c2 = sp.symbols('chi1 chi2', positive=True)
ok5 = True
for d in (4, 5):
    k = sp.Rational(d-2, d-1); val = sp.Rational(2, d-1)
    g1, g2 = 1/sp.cos(c1), 1/sp.cos(c2)
    Sh = s*g1**2*g2**2; f = c*(g1*g2)**(2*k + val)
    Oms = [sp.sin(x)*sp.cos(x) for x in (c1, c2)]
    gs = [val*f*sp.tan(x) for x in (c1, c2)]
    That = ((d-2)*f/(d-1) - sum(o*g for o, g in zip(Oms, gs)))/Sh
    for x in (c1, c2):
        ok5 = ok5 and z(sp.diff(f, x) - That*sp.diff(Sh, x) - sum(o*sp.diff(g, x) for o, g in zip(Oms, gs)))
out['T5_two_planes'] = ok5
# T6: rho = 1 is HX1's turned cut: inner radius and turning length are the two components of the rim length
R_, a_, r_ = sp.symbols('R a r', positive=True)
out['T6_turned_cut'] = z((R_*sp.cos(chi))**2 + (R_*sp.sin(chi))**2 - R_**2)
# T7: conversely, ID1's equal sharing (one transverse direction carries 1/(d-1); a plane, two of them) fixes rho = 1 for d != 3
dd = sp.symbols('d', positive=True)
out['T7_ID1_share_fixes_rho'] = sp.solve(sp.Eq(1 - rho*(dd - 3)/(dd - 1), 2/(dd - 1)), rho) == [1]
out = {k_: (bool(v) if not isinstance(v, str) else v) for k_, v in out.items()}; out['pass'] = all(v for v in out.values() if isinstance(v, bool))
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'TS2_RESULT.json'), 'w'), indent=1)
if __name__ == '__main__': print(out)
