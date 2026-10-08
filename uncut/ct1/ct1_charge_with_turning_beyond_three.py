"""CT1: charge together with turning in d >= 4 directions. A first-law-consistent pure-power candidate exists, and
is refused by the weakly charged solution (Aliev, arXiv:hep-th/0604207, eq. 17); the second-order law that stands.
Exact sympy."""
import json, os, sympy as sp
r, chi, q, s, Qc, th, a_ = sp.symbols('r chi q s Q theta a', positive=True)
out = {}
z = lambda e: sp.simplify(e) == 0
def candidate(d):
    nu = sp.Rational(d-3, d-1)
    gam = 1/sp.cos(chi); c = (d-1)*s/(4*sp.pi); kap = sp.Rational(2, d-1)
    S = s*r**(d-1)*gam**2                                   # TS2 with rho = 1
    M = c*r**(d-2)*gam**2 + c*q**2/r**(d-2)                 # LP2's stored cost added at the inner radius
    J = kap*M*r*sp.tan(chi); Q = q*gam**nu
    Om = sp.sin(chi)*sp.cos(chi)/r
    T, Phi = sp.symbols('T Phi')
    eqs = {x: sp.diff(M, x) - T*sp.diff(S, x) - Om*sp.diff(J, x) - Phi*sp.diff(Q, x) for x in (r, chi, q)}
    sol = sp.solve([eqs[r], eqs[q]], [T, Phi], dict=True)[0]
    factor = sp.simplify(sol[Phi]*Q/(2*c*Q**2/r**(d-2)))    # (energy at second order in the charge)/(static stored cost, same Q)
    return sp.simplify(eqs[chi].subs(sol)), factor
ok_cons = True; fac = {}
for d in (3, 4, 5, 6):
    res, f = candidate(d); ok_cons = ok_cons and z(res); fac[d] = f
out['T1_candidate_is_first_law_consistent'] = ok_cons
# T2: its second-order energy factor: [(d-3) + 2 cos^2]/(d-1) * cos^(2(d-3)/(d-1))
out['T2_candidate_factor'] = all(all(abs(sp.N((fac[d] - ((d-3) + 2*sp.cos(chi)**2)/(d-1)*sp.cos(chi)**sp.Rational(2*(d-3), d-1)).subs(chi, t))) < 1e-12 for t in (0.3, 0.9, 1.3)) for d in fac)
# T3: weakly charged solution, eq. (17): A = -Q/((d-2) r^(d-4) Sigma) (dt - a sin^2 theta dphi); co-turning potential at the horizon
dd = sp.symbols('d', positive=True)
Sigma = r**2 + a_**2*sp.cos(th)**2
At = -Qc/((dd-2)*r**(dd-4)*Sigma); Aphi = -At*a_*sp.sin(th)**2
Omega = a_/(r**2 + a_**2)
PhiH = sp.simplify(-(At + Omega*Aphi))
out['T3_potential_is_uniform_on_the_horizon'] = z(sp.diff(PhiH, th))
out['T3_potential'] = z(PhiH.subs(a_, r*sp.tan(chi)) - Qc*sp.cos(chi)**2/((dd-2)*r**(dd-2)))
true_factor = sp.cos(chi)**2                                  # Phi_H(chi)/Phi_H(0): energy (1/2) Phi Q relative to the static one
# T4: the candidate agrees with it only in three directions
out['T4_agree_d3'] = z(fac[3] - true_factor)
out['T4_refused_d4_d5_d6'] = all(abs(sp.N((fac[d] - true_factor).subs(chi, 1))) > 1e-3 for d in (4, 5, 6))
out['T4_values_at_chi_1'] = str({d: round(float(sp.N(fac[d].subs(chi, 1))), 4) for d in fac}) + ' vs ' + str(round(float(sp.N(true_factor.subs(chi, 1))), 4))
# T5: the weakly charged horizon function (eq. 22) in HL1's form: stored cost enters P as a sum with the turning factor
P = (1 + a_**2/r**2) + q**2*r**(-2*(dd-2))
Delta = r**(dd-2)*(r**2 + a_**2) - sp.Symbol('m')*r**2 + q**2*r**(4-dd)
out['T5_HL1_form'] = z(sp.expand(Delta/r**dd - (P - sp.Symbol('m')*r**(-(dd-2)))))
# T6: the second-order law in one variable u = S / r (count over inner radius):  M = K [ u + (s q)^2 / u ],  K = (d-1)/(4 pi)
ok6 = True
for d in (3, 4, 5, 6):
    gam = 1/sp.cos(chi); c = (d-1)*s/(4*sp.pi); K = sp.Rational(d-1, 4)/sp.pi
    S = s*r**(d-1)*gam**2; u = S/r
    ok6 = ok6 and z(K*u - c*r**(d-2)*gam**2)                                         # neutral turning centre (HD1 b)
    ok6 = ok6 and z(K*(s*q)**2/u - sp.cos(chi)**2*c*q**2/r**(d-2))                   # rim-unit-squared times the static stored cost
    ok6 = ok6 and z((K*(u + (s*q)**2/u)).subs(chi, 0) - (c*r**(d-2) + c*q**2/r**(d-2)))   # static charged centre, exact (LP2)
out['T6_one_variable_form'] = ok6
out = {k_: (bool(v) if not isinstance(v, str) else v) for k_, v in out.items()}; out['pass'] = all(v for v in out.values() if isinstance(v, bool))
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'CT1_RESULT.json'), 'w'), indent=1)
if __name__ == '__main__': print(out)
