"""CF1: why the rim's momentum carries the clock factor. Local quadrature (MT1) + far reading = clock factor x local
reading (AR1-T4, UN2) give E^2 = M0^2 + N^2 p^2. Checked against the published formula, eq. (43) of
Caldarelli, Cognola, Klemm (hep-th/9908022), quoted from the source. Exact sympy."""
import json, os, sympy as sp
m, p, N, R, Q, J, l, h, S = sp.symbols('m_loc p N R Q J l h S', positive=True)
out = {}
z = lambda e: sp.simplify(e) == 0
# T1: local reading in quadrature, then brought to the far reader by the clock factor of the place
E_loc = sp.sqrt(m**2 + p**2)
E_far = N*E_loc
M0 = sp.symbols('M0', positive=True)
out['T1_rule'] = z((E_far**2).subs(m, M0/N) - (M0**2 + N**2*p**2))        # M0 = N m_loc is the far-read rest value
out['T1_flat_is_SP1'] = z((M0**2 + N**2*p**2).subs(N, 1) - (M0**2 + p**2))
# T2: published formula, eq. (43), verbatim in S, J, Q, l
eq43 = S/(4*sp.pi) + sp.pi/(4*S)*(4*J**2 + Q**4) + Q**2/2 + J**2/l**2 + S/(2*sp.pi*l**2)*(Q**2 + S/sp.pi + S**2/(2*sp.pi**2*l**2))
# the line's rules with the ambient clock factor of that space at the rim, N^2 = 1 + R^2/l^2, R = sqrt(S/pi)
Rr = sp.sqrt(S/sp.pi)
M0_line = Rr/2 + Q**2/(2*Rr) + Rr**3/(2*l**2)                             # stored-cost rule (LC1, XU2 with the sign of this space)
line = M0_line**2 + (J/Rr)**2*(1 + Rr**2/l**2)
out['T2_matches_eq43'] = z(sp.expand(line - eq43))
# T3: the ambient clock factor is the neutral space's own: N^2 = 1 + r^2/l^2 is M0's clock with the centre removed
r = sp.symbols('r', positive=True)
Mc = sp.symbols('M', positive=True)
N2_full = 1 - 2*Mc/r + Q**2/r**2 + r**2/l**2
out['T3_ambient_factor'] = z(N2_full.subs({Mc: 0, Q: 0}) - (1 + r**2/l**2)) and z(sp.solve(sp.Eq(N2_full.subs(r, R), 0), Mc)[0] - (R/2 + Q**2/(2*R) + R**3/(2*l**2)))
# T4: expanding space: same statements with 1/l^2 -> -h^2 (continuation used in ET2)
out['T4_continuation'] = z(line.subs(l, sp.I/h) - ((Rr/2 + Q**2/(2*Rr) - h**2*Rr**3/2)**2 + (J/Rr)**2*(1 - h**2*Rr**2)))
# T5: the partner of spin with the factor: Omega = N^2 (J/R^2)/E = N v / R, v = N p / E the rim speed read far
E = sp.sqrt(M0**2 + N**2*(J/R)**2)
out['T5_partner'] = z(sp.diff(E, J) - N*((N*J/R)/E)/R)
out = {k: bool(v) for k, v in out.items()}; out['pass'] = all(out.values())
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'CF1_RESULT.json'), 'w'), indent=1)
if __name__ == '__main__': print(out)
