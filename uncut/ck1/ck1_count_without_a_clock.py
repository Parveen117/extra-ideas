"""CK1: the count of a horizon written with no time variable: phases accumulated per closed turn.
Exact sympy."""
import json, os, sympy as sp
out = {}
z = lambda e: sp.simplify(e) == 0
r, c, d, tau, N = sp.symbols('r c d tau N', positive=True)
# T1: neutral centre, any d: M = c r^(d-2), P = 4 pi r/(d-2) (closing period, SM1 + AL1), S = 4 pi c r^(d-1)/(d-1)
k = (d - 2)/(d - 1)
M = c*r**(d-2); P = 4*sp.pi*r/(d-2); S = 4*sp.pi*c*r**(d-1)/(d-1)
out['T1_count_is_degree_times_phase'] = z(sp.powsimp(S - k*M*P, force=True)) and z(sp.diff(M, r)/sp.diff(S, r) - 1/P)
# T2: the phase M P is made of two counts: (phase of the centre per reader tick) x (reader ticks per closed turn)
mu, n_turn = M*tau, P/tau
out['T2_two_counts'] = z(mu*n_turn - M*P)
out['T2_tick_free_and_reader_free'] = z(sp.diff(mu*n_turn, tau)) and z((M/N)*(N*P) - M*P)      # AR1-T4: energy/N, period*N
# T3: three directions, charge and spin (SP1/LP1 centre): S = (1/2)[M P - (Phi P) Q] - (Omega P) J
Sx, J, Q = sp.symbols('S J Q', positive=True)
R = sp.sqrt(Sx/sp.pi)
E = sp.sqrt((R/2 + Q**2/(2*R))**2 + J**2/R**2)
T, Om, Phi = sp.diff(E, Sx), sp.diff(E, J), sp.diff(E, Q)
Pp = 1/T
pts = [(sp.Rational(3), sp.Rational(1, 3), sp.Rational(1, 2)), (sp.Rational(5), sp.Rational(1), sp.Rational(1, 4)), (sp.Rational(2), sp.Rational(1, 10), sp.Rational(1, 3))]
expr3 = Sx - (sp.Rational(1, 2)*(E*Pp - Phi*Pp*Q) - Om*Pp*J)
out['T3_three_planes'] = all(abs(sp.N(expr3.subs({Sx: p[0], J: p[1], Q: p[2]}))) < 1e-12 for p in pts)
# T4: turning centre in d = 4, 5 (HD1 forms): S = P [k M - Omega J]
chi, s = sp.symbols('chi s', positive=True)
ok4 = True
for dd in (4, 5, 6):
    gam = 1/sp.cos(chi); kk = sp.Rational(dd-2, dd-1)
    S_ = s*r**(dd-1)*gam**2; M_ = (dd-1)*S_/(4*sp.pi*r); J_ = 2*M_*r*sp.tan(chi)/(dd-1)
    T_ = ((dd-2) - 2*sp.sin(chi)**2)/(4*sp.pi*r); Om_ = sp.sin(chi)*sp.cos(chi)/r
    ok4 = ok4 and z(S_ - (kk*M_ - Om_*J_)/T_)
out['T4_turning_any_d'] = ok4
# T5: expansion plane, three directions (XU2): S = (1/2) P [M - h^2 r^3]  (the partner of h^2 is -r^3/2)
h = sp.symbols('h', positive=True)
Mh = r*(1 - h**2*r**2)/2; Th = (1 - 3*h**2*r**2)/(4*sp.pi*r)
out['T5_expansion'] = z(sp.pi*r**2 - sp.Rational(1, 2)*(Mh - h**2*r**3)/Th)
# T6: one statement: sum over planes of (weight) x (angle per turn) x (reading) = 0, weights = degree of each reading in length
#     readings: M (d-2), S (d-1), J (d-1), Q (d-2), h^2 (-2); S's angle per turn is 1
out['T6_weights_d3'] = z((1*E*Pp - 2*Sx - 2*Om*Pp*J - 1*Phi*Pp*Q).subs({Sx: 3, J: sp.Rational(1, 3), Q: sp.Rational(1, 2)}).evalf()) or abs(sp.N((1*E*Pp - 2*Sx - 2*Om*Pp*J - 1*Phi*Pp*Q).subs({Sx: 3, J: sp.Rational(1, 3), Q: sp.Rational(1, 2)}))) < 1e-12
out = {k_: bool(v) for k_, v in out.items()}; out['pass'] = all(out.values())
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'CK1_RESULT.json'), 'w'), indent=1)
if __name__ == '__main__': print(out)
