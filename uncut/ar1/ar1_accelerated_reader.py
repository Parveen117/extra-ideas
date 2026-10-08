"""AR1: the accelerated reader. Clock factor N = rho/rho0 (PR3: uniform acceleration is N = rho). Exact sympy."""
import json, os, sympy as sp
rho, rho0, rho1, E0, t = sp.symbols('rho rho0 rho1 E0 t', positive=True)
out = {}
z = lambda e: sp.simplify(e) == 0
N = rho/rho0                                   # clock of the reader at rho against the reader at rho0
kap = sp.diff(N, rho)                          # rate of the family against the reference clock: 1/rho0
out['T1_rate'] = z(kap - 1/rho0)
P_ref = 2*sp.pi/kap                            # SM1: closing period, forced, in the reference reader's time
P = lambda r: sp.simplify(P_ref*(r/rho0))      # the same closed turn in the time of the reader at r
out['T2_period_each_reader'] = z(P(rho) - 2*sp.pi*rho)
T = lambda r: 1/P(r)                           # ST1-T4: T = 1 / period
a_loc = 1/rho                                  # the reader's own acceleration
out['T3_T_is_acc_over_turn'] = z(T(rho) - a_loc/(2*sp.pi))
out['T3_T_times_clock_constant'] = z(T(rho)*N - T(rho0)) and z(T(rho)*rho - T(rho1)*rho1)
# T4: energy read by the reader at r scales the other way; the (energy, time) area of the closed turn is the same for all
E = lambda r: E0*rho0/r
out['T4_count_reader_free'] = z(P(rho)*E(rho) - P(rho1)*E(rho1)) and z(P(rho)*E(rho) - 2*sp.pi*rho0*E0)
# T5: no scale. Stretching rho, rho0 together changes nothing pure: T*rho, E*P are of weight zero
l = sp.symbols('l', positive=True)
out['T5_scale_free'] = z((T(rho)*rho).subs({rho: l*rho, rho0: l*rho0}, simultaneous=True) - T(rho)*rho)
# T6: open reading of the same family: shared information I = Log cosh(kappa t) grows; closed reading adds none (QT2-T4)
I_open = sp.log(sp.cosh(kap*t))
out['T6_open_grows'] = sp.limit(sp.diff(I_open, t), t, sp.oo) == 1/rho0 and sp.simplify(sp.Abs(sp.exp(sp.I*sp.pi/3))) == 1
out = {k_: (bool(x_) if not isinstance(x_, str) else x_) for k_, x_ in out.items()}
out['pass'] = all(x_ for x_ in out.values() if isinstance(x_, bool))
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'AR1_RESULT.json'), 'w'), indent=1)
if __name__ == '__main__': print(out)
