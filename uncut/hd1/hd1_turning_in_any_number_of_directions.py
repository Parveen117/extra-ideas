"""HD1: turning centres in d = 3..6 directions. The energy-quadrature of SP1 is special to d = 3; the form that holds
in every d; the heat-reading law. Structural inputs are the published family (Emparan & Reall, arXiv:0801.3471,
eqs. 34, 36-38, 45-47), written in the rim reader's coin. Exact sympy."""
import json, os, sympy as sp
out = {}
z = lambda e: sp.simplify(e) == 0
def centre(d, n):
    r, Om = sp.symbols('r Omega', positive=True)
    ch = sp.symbols(f'chi1:{n+1}', positive=True)
    S = (Om/4)*r**(d-1)*sp.prod([1/sp.cos(c)**2 for c in ch])     # (a) each turning plane has rim radius r / cos(chi_i)
    M = (d-1)*S/(4*sp.pi*r)                                        # (b) M = T0(r) S / k, k = (d-2)/(d-1)
    Js = [2*M*r*sp.tan(c)/(d-1) for c in ch]                       # (c) J_i = (2/(d-1)) M R_i v_i
    X = [r] + list(ch)
    T = sp.symbols('T'); Os = sp.symbols(f'O1:{n+1}')
    sol = sp.solve([sp.Eq(sp.diff(M, x), T*sp.diff(S, x) + sum(Os[i]*sp.diff(Js[i], x) for i in range(n))) for x in X], [T] + list(Os), dict=True)[0]
    return r, ch, S, M, Js, sol[T], [sol[o] for o in Os]
law_ok, part_ok = {}, {}
for d, n in ((3, 1), (4, 2), (5, 2), (6, 3)):
    r, ch, S, M, Js, T, Os = centre(d, n)
    v2 = [sp.sin(c)**2 for c in ch]
    law_ok[d] = z(4*sp.pi*r*T - ((d - 2) - 2*sum(v2)))             # heat-reading law
    part_ok[d] = all(z(Os[i]*r - sp.sin(ch[i])*sp.cos(ch[i])) for i in range(n))   # Omega_i = v_i / R_i
out['T1_heat_reading_law_d3_to_d6'] = all(law_ok.values())
out['T2_partner_is_rim_speed_over_rim_radius'] = all(part_ok.values())
# the published temperature of the singly turning member (eq. 38, their d = ours + 1) in the coin
dd, r0, a = sp.symbols('d r0 a', positive=True)
T38 = (2*r0/(r0**2 + a**2) + (dd + 1 - 5)/r0)/(4*sp.pi)
chi = sp.symbols('chi', positive=True)
out['T1_matches_eq38'] = z((4*sp.pi*r0*T38).subs(a, r0*sp.tan(chi)) - ((dd - 2) - 2*sp.sin(chi)**2))
# T3: at fixed count, M = M0 * prod gamma_i^(2/(d-1));  d = 3 is the plain boost of MT1/SP1
for d, n in ((3, 1), (4, 1), (5, 1)):
    r, ch, S, M, Js, T, Os = centre(d, n)
    S0, Om = sp.symbols('S0 Omega', positive=True)
    rsol = sp.solve(sp.Eq(S, S0), r)
    rsol = [s for s in rsol if s.is_positive is not False][0]
    ratio = sp.simplify(M.subs(r, rsol)/M.subs(r, rsol).subs(ch[0], 0))
    out[f'T3_boost_share_d{d}'] = all(abs(sp.N((ratio - sp.cos(ch[0])**sp.Rational(-2, d - 1)).subs(ch[0], t))) < 1e-12 for t in (0.2, 0.7, 1.1))
# T4 (refusal): in four directions the energy-quadrature fails: M^3, not M^2, is linear in J^2 at fixed count
r, ch, S, M, Js, T, Os = centre(4, 1)
S0 = sp.symbols('S0', positive=True)
rs = [s for s in sp.solve(sp.Eq(S, S0), r) if s.is_positive is not False][0]
Mx, Jx = M.subs(r, rs), Js[0].subs(r, rs)
c3 = sp.simplify((Mx**3 - Mx.subs(ch[0], 0)**3)/Jx**2)
c2 = sp.simplify((Mx**2 - Mx.subs(ch[0], 0)**2)/Jx**2)
out['T4_cube_is_linear'] = z(sp.diff(c3, ch[0]))
out['T4_square_is_not'] = not z(sp.diff(c2, ch[0]))
# T5: the surface T = 0 is sum v_i^2 = (d-2)/2. One turning plane reaches it only for d = 3.
out['T5_single_plane_only_d3'] = [d for d in range(3, 9) if sp.Rational(d - 2, 2) < 1] == [3]
# two planes in four directions: v1^2 + v2^2 = 1 ; with r^2 = a b this is the published extremal member
A_, B_ = sp.symbols('a b', positive=True)
out['T5_two_planes_d4'] = z(A_**2/(A_*B_ + A_**2) + B_**2/(A_*B_ + B_**2) - 1)
out = {k: bool(v) for k, v in out.items()}; out['pass'] = all(out.values())
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'HD1_RESULT.json'), 'w'), indent=1)
if __name__ == '__main__': print(out)
