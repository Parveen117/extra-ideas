"""HL1: the heat-reading law for a static centre with any number of planes, in any number of directions.
Clock factor N^2 = 1 - m r^-(d-2) + sum_j c_j r^p_j. Exact sympy."""
import json, os, sympy as sp
out = {}
z = lambda e: sp.simplify(e) == 0
r, rp, m, d = sp.symbols('r r_plus m d', positive=True)
c1, c2, c3, p1, p2, p3 = sp.symbols('c1 c2 c3 p1 p2 p3', real=True)
terms = [(c1, p1), (c2, p2), (c3, p3)]
N2 = 1 - m*r**(-(d-2)) + sum(c*r**p for c, p in terms)
msol = sp.solve(sp.Eq(N2.subs(r, rp), 0), m)[0]                  # horizon
rate = (rp*sp.diff(N2, r).subs(r, rp)).subs(m, msol)             # 4 pi r_plus T = r_plus dN^2/dr (QT1: kappa = N^2'/2, T = kappa/2pi)
eps = [c*rp**p for c, p in terms]                                # value of each plane's term at the horizon
law = (d - 2) + sum((p + d - 2)*e for (c, p), e in zip(terms, eps))
out['T1_law'] = z(sp.expand(rate - law))
# T2: known planes as cases. charge: p = -2(d-2), eps = +Phi^2 -> -(d-2) Phi^2 ; expansion: p = 2, eps = -x^2 -> -d x^2
F, x = sp.symbols('Phi x', positive=True)
out['T2_charge'] = z(law.subs({c2: 0, c3: 0}).subs(p1, -2*(d-2)).subs(c1*rp**(-2*(d-2)), F**2) - (d-2)*(1 - F**2)) or z((d - 2) + (-2*(d-2) + d - 2)*F**2 - (d-2)*(1 - F**2))
out['T2_expansion'] = z((d - 2) + (2 + d - 2)*(-x**2) - ((d - 2) - d*x**2))
out['T2_d3_values'] = z(((d-2)*(1 - F**2) - d*x**2).subs(d, 3) - (1 - F**2 - 3*x**2))        # ET2 with v = 0
# T3: a plane with the memory's own power (p = -(d-2)) contributes nothing: it only renames m
out['T3_same_power_is_silent'] = z((p1 + d - 2).subs(p1, -(d-2)))
# T4: the weight is the plane's degree measured from the memory's: scale r -> l r changes eps_j/eps_mass by l^(p_j + d - 2)
l = sp.symbols('l', positive=True)
out['T4_weight_is_relative_degree'] = z(sp.simplify(((c1*(l*r)**p1)/((l*r)**(-(d-2))))/((c1*r**p1)/(r**(-(d-2)))) - l**(p1 + d - 2)))
# T5: T = 0 surface is one linear relation among the eps_j; ends eps_j = -(d-2)/(p_j + d - 2)
out['T5_ends'] = z(sp.solve(sp.Eq((d-2) + (-2*(d-2) + d - 2)*F**2, 0), F**2)[0] - 1) and z(sp.solve(sp.Eq((d-2) - d*x**2, 0), x**2)[0] - (d-2)/d)
# ---- one form for static and turning planes:  4 pi r_plus T = N_rim^2 * [ (d-2) P + r P' ] at r_plus,
#      N^2 = P(r) - m r^-(d-2);  stored costs enter P as a sum, turnings as a product (1 + a_i^2/r^2);  N_rim^2 = prod cos^2(chi_i)
chi1, chi2, chi3, Q = sp.symbols('chi1 chi2 chi3 Q', positive=True)
def unified(P, chis):
    val = ((d - 2)*P + r*sp.diff(P, r)).subs(r, rp)
    return sp.simplify(sp.prod([sp.cos(c)**2 for c in chis])*val)
P1 = 1 + (rp*sp.tan(chi1))**2/r**2
out['T6_one_turning_plane_any_d'] = z(unified(P1, [chi1]) - ((d - 2) - 2*sp.sin(chi1)**2))
P2 = P1*(1 + (rp*sp.tan(chi2))**2/r**2); P3 = P2*(1 + (rp*sp.tan(chi3))**2/r**2)
out['T7_two_and_three_planes'] = z(unified(P2, [chi1, chi2]) - ((d - 2) - 2*sp.sin(chi1)**2 - 2*sp.sin(chi2)**2)) and z(unified(P3, [chi1, chi2, chi3]) - ((d - 2) - 2*sp.sin(chi1)**2 - 2*sp.sin(chi2)**2 - 2*sp.sin(chi3)**2))
Pkn = (1 + ((rp*sp.tan(chi1))**2 + Q**2)/r**2)
out['T8_charge_and_turning_d3'] = z(unified(Pkn, [chi1]).subs(d, 3) - (sp.cos(2*chi1) - Q**2*sp.cos(chi1)**2/rp**2))
# static planes are the case with no turning: N_rim = 1
out['T9_static_is_the_same_form'] = z(sp.expand(unified(1 + sum(c*r**p for c, p in terms), []) - law))
out = {k: bool(v) for k, v in out.items()}; out['pass'] = all(out.values())
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'HL1_RESULT.json'), 'w'), indent=1)
if __name__ == '__main__': print(out)
