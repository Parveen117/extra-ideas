"""RT1: readers on a rotating system. Exact sympy."""
import json, os, sympy as sp
r, w, P0, E, J, w1, w2, P1, P2, dE, dJ = sp.symbols('r omega P0 E J omega1 omega2 P1 P2 dE dJ', positive=True)
out = {}
z = lambda e: sp.simplify(e) == 0
beta = w*r                                     # rim speed of the co-rotating reader at r
N = sp.sqrt(1 - beta**2)                       # its clock against the reader on the axis
out['T1_one_coin'] = z(N**2 + beta**2 - 1)     # same coin as GR1: N^2 + m = 1 with m = (omega r)^2
# T2: a closed turn shared by the whole system has period P0 on the axis, N P0 for the reader at r
T = 1/(N*P0)
out['T2_T_times_clock_constant'] = z(T*N - 1/P0)
# T3: the unit of the reader at r is N (UN2): 1 on the axis, 0 at omega r = 1: there no count can be read
out['T3_unit_zero_at_rim'] = N.subs(r, 1/w) == 0 and N.subs(r, 0) == 1 and sp.simplify((1/T).subs(r, 1/w)) == 0
# T4: the count is reader-free with the spin plane: energy read at r is (dE - omega dJ)/N
out['T4_count_reader_free'] = z((N*P0)*((dE - w*dJ)/N) - P0*(dE - w*dJ))
# T5: two parts exchanging energy and spin at fixed totals: count stationary for every exchange iff same period and same omega
dS = (P1*(dE - w1*dJ)) - (P2*(dE - w2*dJ))    # part 1 gains (dE, dJ), part 2 loses them (PL1: dS = P dE - (Omega P) dJ)
sol = sp.solve([sp.diff(dS, dE), sp.diff(dS, dJ)], [P1, w1], dict=True)
out['T5_equilibrium_is_rigid'] = sol == [{P1: P2, w1: w2}]
# T6: the same map as the static centre: omega r <-> fall speed; rim <-> horizon (HX1: tanh psi = speed)
psi = sp.atanh(beta)
out['T6_cut_angle'] = z(sp.cosh(psi) - 1/N) and sp.simplify((1/sp.cosh(psi)).subs(r, 1/w)) == 0
out = {k_: (bool(x_) if not isinstance(x_, str) else x_) for k_, x_ in out.items()}
out['pass'] = all(x_ for x_ in out.values() if isinstance(x_, bool))
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'RT1_RESULT.json'), 'w'), indent=1)
if __name__ == '__main__': print(out)
